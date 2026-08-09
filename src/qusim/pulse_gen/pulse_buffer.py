# -*- coding: utf-8 -*-
"""
@author: Meng Wang, Jiheng Duan
"""
import numpy as np
from qutip import *
from qusim.pulse_gen.pulse_config import PulseConfig

def merge_pulse_chan(pulse_buffer_lst: list, pulse: PulseConfig, Hd_i: list):
    """
    Accumulate a drive term ``Hd_i`` into the pulse buffer.

    ``pulse_buffer_lst`` holds three parallel lists: pulse types, qubit indices
    (as strings) and the matching ``[operator, coefficient_array]`` rows. If a row
    already matches both this pulse's type and qubit index, add the new coefficient
    array onto that row; otherwise append a new channel. There is at most one row
    per ``(type, qindex)`` pair, so the first match is the only match.
    """
    types, qindices, rows = pulse_buffer_lst
    qindex = f"{pulse.qindex}"
    for t_index, (chan_type, chan_qindex) in enumerate(zip(types, qindices)):
        if chan_type == pulse.pulse_type and chan_qindex == qindex:
            rows[t_index][1] += Hd_i[1]
            return pulse_buffer_lst

    types.append(pulse.pulse_type)
    qindices.append(qindex)
    rows.append(Hd_i)
    return pulse_buffer_lst



