from qusim.Utils.gate_set import SystemInfo, GateSet
from qusim.PulseGen.pulse_config import PulseConfig, PulseShapeFn
from typing import Optional 
import numpy as np
from copy import deepcopy

class PulseSeqCompiler:

    def __init__(self, gate_sequence: list[list[str]], gate_set: GateSet):
        self.gate_set = gate_set
        self.system_info = self.gate_set.system_info
        self.gate_param = self.gate_set.gate_param
        self.gate_info = self.gate_set.gate_info
        self.gate_seq = gate_sequence
        self.t_moment = self.system_info.t_moment
        
        self.gate_seq_validity()


    def gate_seq_validity(self):
        # The number of gate sequence must equal to the number of qubit
        assert len(self.gate_seq) == self.system_info.num_qubit
        try:
            gseq = np.array(self.gate_seq)
        except:
            raise(ValueError("Gate sequence should be inhomogeneous, each qubit should has the same number of gates"))
        
        pass
    

    def compile_pseq(
            self,
            vritual_Z: Optional[bool] = True,
            init_phase: Optional[np.ndarray[float]] = None
        ):
        self.gate_seq_validity()
        gseq_len = np.shape(self.gate_seq)[1]
        pseq = []; t_delay_buffer = 0
        if init_phase:
            vz_phase_track = init_phase
        else:
            vz_phase_track = np.zeros(self.system_info.num_qubit)
        
        for _i in range(gseq_len):

            num_of_moment_lst = []
            pseq_buffer_lst = []
            for _qidx in range(self.system_info.num_qubit):

                gate_label = self.gate_seq[_qidx][_i]
                
                if gate_label=="I":
                    continue

                if gate_label.find("Z")==0: # add an exam function that exam the qindex
                    """
                    The Z gate label looks like :
                        'Z23.231_0'
                        Z gate on qubit 0, with angle 23.231 degree
                    """
                    slab_index = gate_label.find("_")
                    assert int(gate_label[slab_index+1:])==_qidx
                    if vritual_Z:
                        vz_phase_track[_qidx] -= float(gate_label[1:slab_index-1])
                        continue
                    else:
                        pass

                try:
                    gate_pseq = self.gate_set.dict2pulseconfig(gate_label)
                except:
                    raise(KeyError(f"No such gate in gate_param {gate_label}"))
                try:
                    this_gate_info = self.gate_info[gate_label]
                except:
                    raise(KeyError(f"No such gate in gate_info {gate_label}"))
                
                num_of_moment_lst.append(this_gate_info["num_of_moment"])
                gate_pseq = self.update_t_delay_phase(gate_pseq, t_delay_buffer, vz_phase_track[_qidx])
                pseq_buffer_lst.extend(deepcopy(gate_pseq))

            if len(num_of_moment_lst)>=1:
                if np.max(num_of_moment_lst) == np.min(num_of_moment_lst):
                    pseq.extend(pseq_buffer_lst)
                else:
                    raise(ValueError(f"Imcompatible number of t_moment: gate layer {_i}"))
                t_delay_buffer += np.max(num_of_moment_lst)*self.t_moment
            else:
                pass

        return pseq, t_delay_buffer, vz_phase_track# t_delay_buffer indicates the total gate length


    def update_t_delay_phase(self, gate_pseq: list[PulseConfig], t_delay_new: float, phase_new: float):
        for _i in range(len(gate_pseq)):
            gate_pseq[_i].t_delay += t_delay_new
            if gate_pseq[_i].frequency != 0:
                gate_pseq[_i].phase += phase_new

        return gate_pseq
