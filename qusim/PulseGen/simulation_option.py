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
        init_state: Optional[Qobj] = None,
        dt: Optional[float] = None
    ):
        if dt is not None:
            if dt <= 0:
                raise ValueError("dt must be positive")
            # Prefer dt as source of truth; include endpoint
            simu_point = int(round(simu_time / dt)) + 1
        self.simulation_time = simu_time
        self.simulation_point = simu_point
        self.initial_state = init_state

        self.tlist = np.linspace(0, self.simulation_time, self.simulation_point, dtype=np.float64)

        if self.simulation_point > 1:
            self.dt = self.tlist[1] - self.tlist[0]
            if dt is not None:
                # snap to requested dt to avoid floating drift
                self.dt = dt
        else:
            self.dt = None
