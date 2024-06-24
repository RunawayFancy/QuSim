# -*- coding: utf-8 -*-
"""
@author: Zhen Chen, Jiheng Duan
"""
import os
import pickle
import numpy as np
from sympy import*
from scipy.signal import find_peaks

def idle_w(w1,w2,wcut):     
    ### caculate the detuning off frequency
    x = Symbol('x')
    z = solve((1/(w1-x)+1/(w2-x))/2-1/(w1-wcut))
    z1 = np.abs(z[0]-wcut)
    z2 = np.abs(z[1]-wcut)
    if z1>z2:
        zz = z[1]
    else:
        zz = z[0]
    return zz

def fsim_zhen(u11,u12,u21,u22,u33,result_list, state_001, state_100):
    """
    @author: Zhen Chen
    =====================
    Calculate and return the required parameters of fsim gate

    """
    p0 = state_001.dag()*result_list[0].states[-1].data
    p1 = state_100.dag()*result_list[0].states[-1].data
    print(p1)
    print(p0)
    p0 = np.abs(p0)
    p1 = np.abs(p1)
    theta = np.arctan(p1/p0)[0][0]
    # print(p0)
    # print(p1)
    # print(np.arcsin(p0)[0][0])
    # print(np.arccos(p1)[0][0])
#     phi4 = (u12/(-1j*np.sin(theta)))
    phi2 = (u21/(-1j*np.sin(theta)))
    phi1 = (u22/(np.cos(theta)))
    phi4 = (u11/(np.cos(theta))/phi2)
    phi3 = 1
    phi = u33/phi2/phi4/phi1

    print('theta:', theta)
    print('phi:', phi)
    print('phi1:', phi1)
    print('phi2:', phi2)
    print('phi3:', phi3)
    print('phi4:', phi4)

    return theta, phi, phi1,phi2, phi3, phi4      