# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
from qutip import *

def merge_pulse_chan(pulse_buffer_lst: list, pulse: dict, Hd_i: list):
    index_type = np.where(np.array(pulse_buffer_lst[0])==pulse["type"])[0]
    index_qi = np.where(np.array(pulse_buffer_lst[1])==str(pulse["q_index"]))[0]
    if len(index_type)>0 and len(index_qi)>0:
        t_index = int(np.intersect1d(index_type, index_qi))
        pulse_buffer_lst[2][t_index][1] += Hd_i[1]
    else:
        pulse_buffer_lst[0].append(pulse["type"])
        pulse_buffer_lst[1].append(str(pulse["q_index"]))
        pulse_buffer_lst[2].append(Hd_i)
    return pulse_buffer_lst


