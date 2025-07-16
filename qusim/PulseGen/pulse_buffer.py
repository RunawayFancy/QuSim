# -*- coding: utf-8 -*-
"""
@author: Meng Wang, Jiheng Duan
"""
import numpy as np
from qutip import *
from qusim.PulseGen.pulse_config import PulseConfig

def merge_pulse_chan(pulse_buffer_lst: list, pulse: PulseConfig, Hd_i: list):
    # print(pulse_buffer_lst[0])
    # print(pulse_buffer_lst[1])
    if pulse.pulse_type in pulse_buffer_lst[0]:
        index_type = np.where(np.array(pulse_buffer_lst[0])==pulse.pulse_type)
    else: 
        # print(f'{pulse.pulse_index} miss index type')
        index_type = None
    if f"{pulse.qindex}" in pulse_buffer_lst[1]:
        # print(f"{pulse.qindex} in {pulse_buffer_lst[1]}")
        # index_qi = pulse_buffer_lst[1].index(f"{pulse.qindex}")
        index_qi = np.where(np.array(pulse_buffer_lst[1])==f"{pulse.qindex}")
    else:
        # print(f'{pulse.pulse_index} miss index qi')
        index_qi = None
    # index_type = np.where(np.array(pulse_buffer_lst[0])==pulse.pulse_type)[0]
    # index_qi = np.where(np.array(pulse_buffer_lst[1])==str(pulse.qindex))[0]
    if index_type and index_qi:
        # print([index_type, index_qi])
        # print(np.intersect1d(index_type, index_qi))
        t_index = int(np.intersect1d(index_type, index_qi))
        pulse_buffer_lst[2][t_index][1] += Hd_i[1]
    else:
        pulse_buffer_lst[0].append(pulse.pulse_type)
        pulse_buffer_lst[1].append(f"{pulse.qindex}")
        pulse_buffer_lst[2].append(Hd_i)
    # print('---------------')
    return pulse_buffer_lst



