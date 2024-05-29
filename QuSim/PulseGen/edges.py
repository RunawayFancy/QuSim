# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np 

def max_min(A,B):
    return max(A,B), min(A,B)

def rasing_t(tlist, t_delay, _t_width, a, b):
    t_width = max(1e-1, _t_width)
    A = 2*np.pi/t_width *(b-a)
    B = 2*np.pi/2 *(a-1) - A * t_delay
    return A * tlist + B

def lowering_t(tlist, t_delay, _t_width, t_plateau, a, b):
    t_width = max(1e-1, _t_width)
    A = - 2*np.pi/t_width *(b-a)
    B = 2*np.pi/2 *(a-1) - A * (t_delay + t_width + t_plateau)
    return A * tlist + B

# Revise
def cosine_edge(tlist, t_delay, _t_width, t_plateau, _a = 0, _b = 1):
    t_width = max(1e-1, _t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    b, a = max_min(_a, _b)

    Y += ( tlist < t_delay ) * 0
    Y += ( (tlist >= t_delay) & (tlist < t_delay + t_width/2) ) * (1 - np.cos(np.pi * (tlist - t_delay) /(t_width/2)))
    # (np.cos(rasing_t(tlist, t_delay, t_width, a, b)) + 1)
    Y += ((tlist >= t_delay + t_width/2) & (tlist <= t_delay + t_width/2 + t_plateau)) * 2
    Y += ((tlist > t_delay + t_width/2 + t_plateau) & (tlist <= t_delay + t_width + t_plateau)) * (1 - np.cos(np.pi * (tlist - t_delay - t_plateau) / (t_width/2)))
    # (np.cos(lowering_t(tlist, t_delay, t_width, t_plateau, a, b)) + 1)
    Y += ( tlist > (t_delay + t_width + t_plateau) ) * 0
    return Y/np.max(Y)

def hcosh_edge(tlist, t_delay, _t_width, t_plateau):
    zero1, zero2 = np.log(2-np.sqrt(3)), np.log(2+np.sqrt(3))
    t_width = max(1e-1, _t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)

    x1 = zero1 * ((tlist - t_delay)/(0.5 * t_width) - 1)
    x2= zero2 * ((tlist - t_delay - t_width * 0.5 - t_plateau)/(0.5 * t_width))

    Y += ( tlist < t_delay ) * 0
    Y += ( (tlist >= t_delay) & (tlist <= t_delay + t_width/2) ) * ( - np.cosh(x1) + 2)
    Y += ((tlist > t_delay + t_width/2) & (tlist < t_delay + t_width/2 + t_plateau)) * 1
    Y += ((tlist >= t_delay + t_width/2 + t_plateau) & (tlist <= t_delay + t_width + t_plateau)) * ( - np.cosh(x2) + 2)
    Y += ( tlist > (t_delay + t_width + t_plateau) ) * 0
    return Y

def tanh_edge(tlist, t_delay, _t_width, t_plateau, epsilon):
    t_width = max(1e-1, _t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    d = t_delay
    w = t_width
    p = t_plateau
    tau = w/2

    Y += ( tlist < d ) * 0
    Y += ( (( tlist >= d ) & ( tlist <= d + w/2)  )+ 0j ) * np.tanh(2*epsilon*(tlist-d)/tau)
    Y += ((tlist > d + w/2) & (tlist < d + w/2 + p)) * 1
    Y += ( (( tlist >= d + w/2 + p ) & ( tlist <= d + w + p )  )+ 0j) * np.tanh(2*epsilon*(2*tau + p -(tlist-d))/tau)
    Y += ( tlist > (d + w + p) ) * 0
    return Y

def hyper_edge(tlist, t_delay, _t_width, t_plateau, epsilon):
    t_width = max(1e-1, _t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    w = t_width
    d = t_delay
    p = t_plateau
    tau = w

    Y += ( tlist < d ) * 0
    Y += (np.cosh(epsilon/2) - np.cosh(epsilon*((tlist-d)/tau-1/2)))*( (( (tlist-d) >= 1e-11 ) & ( (tlist-d-w/2) <= 1e-11 ) )+ 0j)
    Y += (np.cosh(epsilon/2) - np.cosh(epsilon*((tlist-d-p)/tau-1/2)))*( (( (tlist - d - w/2 - p) >= 1e-11 ) & ( (tlist - d - w - p) <= 1e-11 ) )+ 0j)
    Y/=np.max(Y)
    Y += ((tlist > d + w/2) & (tlist < d + w/2 + p)) * 1
    Y += ( tlist > (d + w + p) ) * 0
    return Y