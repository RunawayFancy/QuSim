# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""

from qusim.PulseGen.simulation_option import SimulationOption 
import numpy as np
import numpy.fft as fft
import inspect
from typing import Literal, Optional, Iterable, Callable

class NoiseTimeConfig:
    """
    simopt: SimulationOption
    tseg: Optional[float]
    tstart: Optional[float]
    tstop: Optional[float]

    time_dependent: bool = False,
        True: The noise series will become time dependent.
            It requires tranfofn.

        False: ~ will be time-independent
            It does not require tranfofn.

    tranfofn: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None
        Transformation function. When the noise is time-dependent, the transformation fuction will use `tlist` and `waveform` as input to generate a time dependent noise sequence.

        tranfofn input
            tranfofn(tlist: np.ndarray, waveform, np.ndarray) -> np.ndarray
            
            assert tlist.shape ==waveform.shape
    """
    def __init__(self,
        simopt: SimulationOption,
        tseg: Optional[float] = None,
        tstart: Optional[float] = None,
        tstop: Optional[float] = None,
        time_dependent: bool = False,
        tranfofn: Optional[Callable] = None,
        tdbasefn: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None
    ):
        self.simopt = simopt
        if tstop:
            self.tstop = tstop
        else:
            self.tstop = self.simopt.simulation_time

        if tstart:
            self.tstart = tstart
        else:
            self.tstart = 0
        
        if tseg:
            self.tseg = tseg
        else:
            self.tseg = self.tstop
        self.time_dependent = time_dependent
        
        self.tranfofn = tranfofn

        if tdbasefn is not None:
            # Ensure tranfofn is callable
            if not callable(tdbasefn):
                raise ValueError("tranfofn should be a callable function or None.")
            
            # Check the signature of tranfofn
            signature = inspect.signature(tdbasefn)
            parameters = list(signature.parameters.keys())
            if len(parameters) != 2:
                raise ValueError("tdbasefn should accept exactly two parameters: waveform and tlist.")
        
        self.tdbasefn = tdbasefn


class GaussianNoiseConfig:
    """
    mean: float
    std: float
    noise_time_config: NoiseTimeConfig
    """
    def __init__(self,
        noise_time_config: NoiseTimeConfig,
        mean: float = 0,
        std: float = 0.1,
        methods: Literal['sum', 'multiply'] = 'sum'
    ):
        self.noise_time_config = noise_time_config
        self.mean = mean
        self.std = std
        self.type = 'ga'
        self.methods = methods
    
    def trigger(self, seg_length: int) -> np.ndarray:
        return np.random.normal(self.mean, self.std, seg_length)


class RandomTeleNoiseConfig:
    """
    noise_time_config: NoiseTimeConfig
    high_val: float
    low_val: float
    switch_prob: float
    """
    def __init__(self,
        noise_time_config: NoiseTimeConfig,
        high_val: float = 1,
        low_val: float = 0,
        switch_prob: float = 0.5,
        methods: Literal['sum', 'multiply'] = 'sum'
    ):
        self.noise_time_config = noise_time_config
        self.high_val = high_val
        self.low_val = low_val
        self.switch_prob = switch_prob
        self.type = 'rt'
        self.methods = methods

    def trigger(self, seg_length: int) -> np.ndarray:
        return np.random.choice([self.low_val, self.high_val], seg_length, p=[self.switch_prob, 1-self.switch_prob])


class JNNoiseConfig:
    """    
    noise_time_config: NoiseTimeConfig
    temperature: float = 295 # Unit in Kelvin
    resistance: float = 1 # Unit in Ohm
    """
    def __init__(self,
        noise_time_config: NoiseTimeConfig,
        temperature: float = 295, # Unit in Kelvin
        resistance: float = 1, # Unit in Ohm
        methods: Literal['sum', 'multiply'] = 'sum'             
    ):
        self.noise_time_config = noise_time_config
        self.temperature = temperature
        self.resistance = resistance
        self.kb = __KB__
        self.type = 'jn'
        self.methods = methods

    def trigger(self, seg_length: int) -> np.ndarray:
        V_noise_rms = np.sqrt(4 * self.kb * self.temperature * self.resistance)
        return np.random.normal(0, V_noise_rms, seg_length)
    

class OneOverFNoiseConfig:
    """
    noise_time_config: NoiseTimeConfig
    alpha: float = 0.9 # 1/f^alpha
    scale: float = 1 # scale/f^alpha
    """
    def __init__(self,
        noise_time_config: NoiseTimeConfig,
        hfreq: float = 1e10,
        lfreq: float = 1e-2,
        alpha: float = 0.9, # 1/f^alpha
        amp: float = 1.0,
        scale: float = 1e10, # scales the ifft noise time series
        normalized: Literal['no', 'std', 'max'] = 'no', # Normalized the ifft signal with std
        methods: Literal['sum', 'multiply'] = 'sum'
    ):
        self.noise_time_config = noise_time_config
        self.alpha = alpha
        self.amp = amp
        self.scale = scale
        self.cutoff_freq_high = hfreq
        self.cutoff_freq_low = lfreq
        self.normalized = normalized
        self.type = '1/f'
        self.methods = methods
    
    def trigger(self, N: int) -> np.ndarray:
        """
        Unfinished.
        """
        dt = self.noise_time_config.simopt.dt*1e-9
        T =  self.noise_time_config.simopt.simulation_time*1e-9
        freqs = np.fft.fftfreq(N, d=dt)
        freq_band_mask = (np.abs(freqs) >= self.cutoff_freq_low) & (np.abs(freqs) <= self.cutoff_freq_high)
        S_f = np.zeros(N)
        mm = self.scale * self.amp
        S_f[freq_band_mask] = mm * np.abs(freqs[freq_band_mask]) ** (-self.alpha)
        phases = np.random.uniform(0, 2 * np.pi, N)
        A_f = S_f * np.exp(1j * phases)
        A_f = np.fft.ifftshift(A_f)
        A_f = np.fft.fftshift(A_f)
        noise_arr = np.fft.ifft(A_f).real
        
        return noise_arr

__KB__= 1.38e-23  # Boltzmann constant in J/K