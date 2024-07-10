from qutip import Qobj
from typing import Optional
import numpy as np

class SimulationOption():
    simulation_time: float
    simulation_point: int
    initial_state: Optional[Qobj] = None
    tlist: np.ndarray

    def __init__(self,
        simu_time: float,
        simu_point: int,
        init_state: Optional[Qobj] = None
    ):
        self.simulation_time = simu_time
        self.simulation_point = simu_point
        self.initial_state = init_state
        
        self.tlist = np.linspace(0, self.simulation_time, self.simulation_point)