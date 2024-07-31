from qusim.Utils.gate_set import SystemInfo, GateSet
from typing import Optional 
import numpy as np


PI = np.pi


class PulseSeqCompiler:

    def __init__(self, gate_sequence: list[list[str]], gate_set: GateSet):
        self.gate_set = gate_set
        self.system_info = self.gate_set.system_info
        self.gate_param = self.gate_set.gate_param
        self.gate_seq = gate_sequence
        self.t_moment = self.system_info.t_moment


        self.gate_seq_validity()
        
        pass

    def gate_seq_validity(self):
        # The number of gate sequence must equal to the number of qubit
        assert len(self.gate_seq) == self.system_info.num_qubit
        


        pass
    
    def compile_pseq(
            self,
            vritual_Z: Optional[bool] = True
        ):
        """
        The self.gate_sequence looks like
        [
        ["X90_0", "I", ]
        ]
        
        """
        
        
        pass