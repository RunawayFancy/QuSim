# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np


PI = np.pi
Jres = 0.1e-3

def Jexchange(v_barrier: np.ndarray, lever_arm: float = 1, J_res: float = Jres):
    return J_res*np.exp(2*lever_arm*v_barrier) # unit in GHz


def J2Vbarrier(J_exchange: np.ndarray, lever_arm: float = 1, J_res: float = Jres):
    if np.min(J_exchange/J_res + 1)<=0:
        print((J_exchange/J_res + 1)<0)
    return 1/(2*lever_arm) * np.log(J_exchange/J_res + 1)
        

def dJdV(v_barrier: np.ndarray, lever_arm = 1, J_res = Jres):
    return J_res*2*lever_arm*np.exp(2*lever_arm*v_barrier) # unit in GHz


def tdbase_charge_noise(tlist: np.ndarray, waveform: np.ndarray):
    Vbarrier_arr = J2Vbarrier(waveform+Jres)
    sensitivity_arr = dJdV(Vbarrier_arr)
    return sensitivity_arr


def tranfofn_charge_noise(v_barrier: np.ndarray, lever_arm: float = 1):
    return np.exp(2*lever_arm*v_barrier)