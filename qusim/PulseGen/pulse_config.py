# -*- coding: utf-8 -*-
"""
@author: Meng Wang, Jiheng Duan
"""

import qusim.PulseGen.edges as edges
# from qusim.PulseGen.pulse_shape import PulseShapeFn
from qusim.PulseGen.noise_gen import noise_gen
from qusim.Instruments.tools import grad
from qusim.PulseGen.simulation_option import SimulationOption 
from qusim.PulseGen.noise_config import * 

from collections import namedtuple
from collections import defaultdict as ddict
from typing import Literal, Optional, Iterable
from enum import Enum

from copy import deepcopy
import numpy as np
from numpy import pi as PI

import socket
import pickle


def cosine(tlist: np.ndarray, pulse: 'PulseConfig'):
    # If the width of the pulse is small, then return a square pulse.
    if np.abs(pulse.t_width) < 1e-5: return square(tlist, pulse)
    return pulse.amplitude * edges.cosine_edge(tlist, pulse.t_delay, pulse.t_width, pulse.t_plateau)

def hcosine(tlist: np.ndarray, pulse: 'PulseConfig'): # Half cosine
    if np.abs(pulse.t_width) < 1e-5: return square(tlist, pulse)
    a = 0.3; b = 1
    return pulse.amplitude * edges.cosine_edge(tlist, pulse.t_delay, pulse.t_width, pulse.t_plateau, a, b)

def cosh(tlist: np.ndarray, pulse: 'PulseConfig'):
    if np.abs(pulse.t_width) < 1e-5: return square(tlist, pulse)
    return pulse.amplitude * edges.hcosh_edge(tlist, pulse.t_delay, pulse.t_width, pulse.t_plateau)

def square(tlist: np.ndarray, pulse: 'PulseConfig'):
    return pulse.amplitude * ( ((tlist >= pulse.t_delay) & (tlist <= pulse.t_delay + pulse.t_width + pulse.t_plateau) ) + 0j )

def tanh(tlist: np.ndarray, pulse: 'PulseConfig'):
    return pulse.amplitude * edges.tanh_edge(tlist, pulse.t_delay, pulse.t_width, pulse.t_plateau, pulse.epsilon)

def hyper(tlist: np.ndarray, pulse: 'PulseConfig'):
    return pulse.amplitude * edges.hyper_edge(tlist, pulse.t_delay, pulse.t_width, pulse.t_plateau, pulse.epsilon)

# class Cosine(Enum):
#     FN = cosine
#     LABEL = "cosine"

# class HCosine(Enum):
#     FN = hcosine
#     LABEL = "hcosine"

# class Cosh(Enum):
#     FN = cosh
#     LABEL = "cosh"

# class Square(Enum):
#     FN = square
#     LABEL = "square"

# class Tanh(Enum):
#     FN = tanh
#     LABEL = "tanh"

# class Hyper(Enum):
#     FN = hyper
#     LABEL = "hyper"

class PulseShapeFn(Enum):
    """
    Pulse Shape Function Enum Class
    
    Usage: 
        - `PulseShapeFn.SQUARE` <=> `square`
        - `PulseShapeFn.SQUARE()` <=> `square()`
        - pass in `PulseShapeFn.SQUARE` as a parameter, e.g. `pulse_shape = PulseShapeFn.SQUARE`
    """
    SQUARE = square
    COSINE = cosine
    HCOSINE = hcosine
    COSH = cosh
    TANH = tanh
    HYPER = hyper

    # SQUARE = Square
    # COSINE = Cosine
    # HCOSINE = HCosine
    # COSH = Cosh
    # TANH = Tanh
    # HYPER = Hyper

# noise_chan1 = [
#     {
#     'type': 'white',
#     'switch': 'off',
#     'std': 0.3,
#     }
# ]

# drag = DRAGConfig(1, 0) when initialized
# drag.scale, drag[0] same thing
# drag.delta, drag[1] same thing
DRAGConfig = namedtuple('DRAGConfig', ['scale', 'delta'], defaults=[0.0, 1.0])
DRAGConfig.__doc__ = """0 <= scale <= 1, delta != 0"""

class PulseConfig():
    """Single Pulse Configuration

    Args:
        pulse_index (int): _description_
        pulse_type (Literal[&#39;XY&#39;, &#39;Z&#39;, &#39;INT&#39;]): _description_
        pulse_shape (PulseShape): e.g. PulseShape.SQUARE
        t_delay (float): _description_
        t_width (float): _description_
        t_plateau (float): _description_
        qindex (int | list): _description_
        phase (float, optional): _description_. Defaults to 0.
        frequency (float, optional): _description_. Defaults to 0.
        amplitude (float, optional): _description_. Defaults to 0.
        noise (Optional[dict], optional): _description_. Defaults to None.
        epsilon (float, optional): _description_. Defaults to 1.
        frequency_detuning (float, optional): _description_. Defaults to 0.
        DRAG_scale (float, optional): _description_. Defaults to 0.
        DRAG_delta (float, optional): _description_. Defaults to 0.
    """

    pulse_index: int
    pulse_type: Literal['XY', 'Z', 'INT']
    pulse_shape: PulseShapeFn

    t_delay: float
    t_width: float
    t_plateau: float
    qindex: int|list

    phase: float = 0
    frequency: float = 0
    amplitude: float = 0
    offset: Optional[float] = None
    noise: Optional[list] = None
    epsilon: float = 1
    frequency_detuning: Optional[float] = None
    
    DRAG_config_list: Optional[Iterable[DRAGConfig]] = None

    predistortion: Optional[list] = None

    def __init__(self,
        pulse_index: int,
        pulse_type: Literal['XY', 'Z', 'INT'],
        pulse_shape: PulseShapeFn,
    
        t_delay: float,
        t_width: float,
        t_plateau: float,
        qindex: int|list,

        phase: float = 0,
        frequency: float = 0,
        amplitude: float = 0,
        offset: Optional[float] = None,
        noise: Optional[list[GaussianNoiseConfig|RandomTeleNoiseConfig|JNNoiseConfig|OneOverFNoiseConfig]] = None,
        epsilon: float = 1,
        frequency_detuning: Optional[float] = None,
        
        DRAG_config_list: Optional[Iterable[DRAGConfig]] = None,

        predistortion: Optional[list] = None
    ):
        self.qindex = qindex
        
        self.pulse_index = pulse_index
        self.pulse_type = pulse_type
        self.pulse_shape = pulse_shape
        self.t_delay = t_delay
        self.t_width = t_width
        self.t_plateau = t_plateau

        self.phase = phase
        self.frequency = frequency
        self.amplitude = amplitude
        self.offset = offset
        self.noise = noise
        self.epsilon = epsilon

        self.frequency_detuning = frequency_detuning
        self.DRAG_config_list = DRAG_config_list

        self.predistortion = predistortion
    
    def __str__(self):
        return \
            f'name={self.pulse_type}{self.qindex},' \
            f'pid={self.pulse_index},'              \
            f'shp={self.pulse_shape.__name__}'
        
    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other: 'PulseConfig'):
        return self.qindex == other.qindex and self.pulse_type == other.pulse_type

    def __lt__(self, other: 'PulseConfig'):
        if self.qindex == other.qindex:
            return self.pulse_type < other.pulse_type

        return self.qindex < other.qindex

    def send2plot(self, host: str = '127.0.0.1', port: int = 63243) -> None:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
        s.send(pickle.dumps(self))
        s.close()

    def export2dict(self) -> dict:
        pulse_param_dict = dict()
        pulse_param_dict["pulse_index"] = self.pulse_index
        pulse_param_dict["pulse_type"] = self.pulse_type
        pulse_param_dict["pulse_shape"] = self.pulse_shape.__name__
        pulse_param_dict["t_delay"] = self.t_delay
        pulse_param_dict["t_width"] = self.t_width
        pulse_param_dict["t_plateau"] = self.t_plateau
        pulse_param_dict["qindex"] = self.qindex
        pulse_param_dict["phase"] = self.phase
        pulse_param_dict["frequency"] = self.frequency
        pulse_param_dict["amplitude"] = self.amplitude
        pulse_param_dict["offset"] = self.offset
        pulse_param_dict["noise"] = self.noise
        pulse_param_dict["epsilon"] = self.epsilon
        pulse_param_dict["frequency_detuning"] = self.frequency_detuning
        pulse_param_dict["DRAG_config_list"] = self.DRAG_config_list
        pulse_param_dict["predistortion"]= self.predistortion

        return pulse_param_dict
        
    def get_pulse(self, sim_opts: SimulationOption) -> np.ndarray:
        """
        Get pulse (1D np.array)

        pulse = Re [waveform * carrier]
        """
        self.amplitude*=2*PI
        tlist = sim_opts.tlist
        delta_t = sim_opts.simulation_time / sim_opts.simulation_point
        
        # Add envolope 
        drive_pulse = self.pulse_shape(tlist, self)

        # Add DRAG correction
        # Multiple DRAG, define a new function
        if self.DRAG_config_list:
            
            num_drag = len(self.DRAG_config_list)
            for i in range(num_drag):
                drive_pulse += self.DRAG(drive_pulse, self.DRAG_config_list[i], delta_t)
        
        # Add frequency deturning
        if self.frequency_detuning:
            t_center = (self.t_width + self.t_plateau)/2 + self.t_delay
            drive_pulse *= np.exp(1j*2*PI * self.frequency_detuning * (tlist - t_center))
        
        # Add carrier
        drive_pulse = np.real(drive_pulse * self.carrier(tlist))
        
        # Add noise
        if self.noise:
            drive_pulse = self.add_noise(drive_pulse, sim_opts)

        # Add offset
        if self.offset:
            drive_pulse += self.add_offset(tlist)
        self.amplitude/=2*PI
        return drive_pulse

    def carrier(self, tlist: np.ndarray) -> np.ndarray:
        return np.exp(-1j * (2*PI * self.frequency * tlist + self.phase))
    
    def DRAG(self, ylist: np.ndarray, DRAG_config: DRAGConfig, delta_t: float) -> np.ndarray:
        assert 0 <= DRAG_config.scale <= 1, "DRAG scale must be in [0, 1]"
        assert DRAG_config.delta != 0, "DRAG delta can not be zero"
        assert delta_t > 0, "Delta t can not be equal and small than zero"
        
        return -1j*grad(ylist)/delta_t * DRAG_config.scale/(2*PI * DRAG_config.delta)
    
    def add_noise(self, ylist: np.ndarray, sim_opts: SimulationOption):
        for noise_config in self.noise:
            noise = noise_gen(noise_config, deepcopy(ylist))
            ylist += np.real(noise)

        return ylist

    def add_offset(self, tlist: np.ndarray) -> np.ndarray:
        pulse = deepcopy(self)
        pulse.amplitude = pulse.offset * 2*PI
        offset_pulse = PulseShapeFn.SQUARE(tlist, pulse)

        return np.real(offset_pulse)
        

# used for testing other functions
# var name begins with __ so that 
# it is not imported when `import *`
__TEST_PULSE__ = PulseConfig(
    # TODO: add the right params
    pulse_index = 0,
    pulse_type = 'INT',
    pulse_shape = PulseShapeFn.COSINE,
    t_delay=0,
    t_width=4,
    t_plateau=0,
    frequency=3,
    phase=0,
    amplitude=0.5,
    qindex=0
)
__TEST_PULSE_2__ = PulseConfig(
    pulse_index=1, 
    pulse_shape=PulseShapeFn.COSINE,
    pulse_type="XY",
    t_delay= 0,
    t_width= 10,
    t_plateau=0,
    frequency= 6.3,
    amplitude= 0.1,
    qindex= 1,
    DRAG_config_list=[DRAGConfig(scale=1, delta= -0.19)]
)
__TEST_PULSE_3__ = PulseConfig(
    pulse_index=1, 
    pulse_shape=PulseShapeFn.COSINE,
    pulse_type="XY",
    t_delay= 0,
    t_width= 10,
    t_plateau=0,
    frequency= 6.3,
    amplitude= 0.1,
    qindex= 1,
    offset= 0.2
)
__TEST_SIM_OPT2__ = SimulationOption(
    simu_point= 100000,
    simu_time= 10,
)

if __name__ == '__main__':
    print(__TEST_PULSE__)
    