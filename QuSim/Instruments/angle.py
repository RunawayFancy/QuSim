# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
from qutip import isket


def cal_angle(result_list, simulation_option):
    angle_list = []
    for ii, res in enumerate(result_list):
        angle = np.angle(simulation_option["initial_state"][ii].dag() * res.states[-1])[0][0]
        angle_list.append(angle)
    return (angle_list[3] - angle_list[2] - angle_list[1] - angle_list[0])

        
def get_angle(tstate, result):
    spt = result.states[-1]
    if isket(spt):
        if isket(tstate):
            angle = np.angle(tstate.dag() * spt)
    # else:
    #     if isket(tstate):
    #         length = spt.dims[0][0]
    #         # print(length)
    #         # print(tstate.dims)
    #         # print(Qobj(np.insert(np.angle(spt.full()[0][1:length]), 0 ,0)).dims)
    #         angle_qobj = Qobj(np.insert(np.angle(spt.full()[0][1:length]), 0 ,0)).dag() * tstate
    #         angle = angle_qobj.full()[0][0]
    #     else: angle = None
    else: angle = None
    # Normalied the angle to [0, 2pi] range
    if angle is not None:
        angle = angle % (2 * np.pi)
    return angle