# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""

from qusim.PulseGen.simulation_option import SimulationOption 
import numpy as np
import numpy.fft as fft
import inspect
from typing import Literal, Optional, Iterable, Callable
from scipy.signal import welch
import matplotlib.pyplot as plt

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

            self.time_dependent = True
            
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
        amp: float = 1,
        methods: Literal['sum', 'multiply'] = 'sum'
    ):
        self.noise_time_config = noise_time_config
        self.mean = mean
        self.std = std
        self.amp = amp
        self.type = 'ga'
        self.methods = methods
    

    def trigger(self, seg_length: int) -> np.ndarray:
        return self.amp*np.random.normal(self.mean, self.std, seg_length)


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
        lfreq: float = 1, # Hz
        alpha: float = 0.9, # 1/f^alpha
        scale: float = 1, # scales the ifft noise time series
        methods: Literal['sum', 'multiply'] = 'sum',
        # normalized: Literal['no', 'std', 'max'] = 'no', # Normalized the ifft signal with std
        
    ):
        self.noise_time_config = noise_time_config
        self.npts = self.noise_time_config.simopt.simulation_point
        self.t_total = self.noise_time_config.simopt.simulation_time * 1e-9 # unit in seconds
        self.alpha = alpha
        self.scale = scale

        self.cutoff_freq_high = self.npts/2/self.t_total
        self.cutoff_freq_low = lfreq
        
        self.ml = self.cutoff_freq_low * self.t_total + 1
        m_series = np.linspace(1, self.npts, self.npts)
        m_series[0] = self.ml
        self.m_series = m_series

        self.sampling_freq = self.npts/self.t_total

        self.methods = methods
        # self.normalized = normalized

        self.type = '1/f'
    

    def trigger(self, N: int) -> np.ndarray:
        assert N == self.npts

        # if self.skip_freq_mask:
        #     S_f = self.scale * np.abs(self.freqs) ** (-self.alpha)
        # else:
        #     freq_band_mask = (np.abs(self.freqs) >= self.cutoff_freq_low) & (np.abs(self.freqs) <= self.cutoff_freq_high)
        #     S_f = np.zeros(N)
        #     S_f[freq_band_mask] = self.scale * np.abs(self.freqs[freq_band_mask]) ** (-self.alpha)

        # Calculate the noise fourier coefficient

        S_f = self.scale * np.abs((self.m_series - 1)/self.t_total) ** (-self.alpha) * N / self.t_total

        # phases = np.random.uniform(0, 2 * np.pi, N)
        gaussian_white_noise_seq = np.random.normal(0,1,N)
        fft_gaussian_white_noise_seq = np.fft.fft(gaussian_white_noise_seq)
        A_f = np.sqrt(S_f) * fft_gaussian_white_noise_seq
        A_f = np.fft.ifftshift(A_f)
        A_f = np.fft.fftshift(A_f)
        noise_arr = np.fft.ifft(A_f).real
        
        return noise_arr
    

    def show(self):
        noise = self.trigger(self.npts)
        frequencies, psd = welch(noise, fs=self.sampling_freq, noverlap=0, nperseg=10000)
        num_bins=100

        fig, axes = plt.subplots(3, 1, figsize=(9, 9))

        # plt.plot(noise_old, label=f'1/f^{alpha} noise')
        axes[0].plot(self.noise_time_config.simopt.tlist*1e-9, noise, label=f'1/f^{self.alpha} noise')
        axes[0].set_title('Time Series with 1/f^alpha Power Spectrum')
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel(r'$\delta V_B$ (mV)')
        axes[0].legend()

        # Plotting the PSD
        axes[1].loglog(frequencies, psd)
        axes[1].set_title('PSD of the Noise Series')
        axes[1].set_xlabel('Frequency (Hz)')
        axes[1].set_ylabel(r'PSD ($mV^2$/Hz)')
        axes[1].grid(True)

        axes[2].hist(noise, bins=num_bins, edgecolor='black', alpha=0.7)
        axes[2].set_title(f'Histogram of the noise time series, num_bins={num_bins}')
        axes[2].set_xlabel(r'$\delta V_B$ (mV)')
        axes[2].set_ylabel('Counts')
        axes[2].grid(True)

        plt.tight_layout()

        pass
    
__KB__= 1.38e-23  # Boltzmann constant in J/K