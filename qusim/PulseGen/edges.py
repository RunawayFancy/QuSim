# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np 
from typing import Literal

def max_min(A,B):

    return max(A,B), min(A,B)


def rasing_t(tlist: np.ndarray, 
            t_delay: float,
            t_width: float,
            a: float,
            b: float):
    
    t_width = max(1e-1, t_width)
    A = 2*np.pi/t_width *(b-a)
    B = 2*np.pi/2 *(a-1) - A * t_delay

    return A * tlist + B


def lowering_t(tlist: np.ndarray, 
               t_delay: float, 
               t_width: float, 
               t_plateau: float, 
               a: float, 
               b: float):
    
    t_width = max(1e-1, t_width)
    A = - 2*np.pi/t_width *(b-a)
    B = 2*np.pi/2 *(a-1) - A * (t_delay + t_width + t_plateau)

    return A * tlist + B


def square_edge(tlist: np.ndarray,
                t_delay: float,
                t_plateau: float):
    d = t_delay
    p = t_plateau

    return ( ((tlist >= d) & (tlist <= d + p) ))


# Revise
def cosine_edge(tlist: np.ndarray,
                t_delay: float,
                t_width: float,
                t_plateau: float,
                _a: float = 0,
                _b:float = 1):
    t_width = max(1e-1, t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    b, a = max_min(_a, _b)
    w = t_width
    d = t_delay
    p = t_plateau

    # print(np.shape(tlist))
    # print(np.shape(Y))
    # print(np.shape(( tlist < t_delay ) * 0))

    Y += ( tlist < d ) * 0
    Y += ( (tlist >= d) & (tlist < d + w/2) ) * (1 - np.cos(np.pi * (tlist - d) /(w/2)))
    # (np.cos(rasing_t(tlist, t_delay, t_width, a, b)) + 1)
    Y += ((tlist >= d + w/2) & (tlist <= d + w/2 + p)) * 2
    Y += ((tlist > d + w/2 + p) & (tlist <= d + w + p)) * (1 - np.cos(np.pi * (tlist - d - p) / (w/2)))
    # (np.cos(lowering_t(tlist, t_delay, t_width, t_plateau, a, b)) + 1)
    Y += ( tlist > (d + w + p) ) * 0

    return Y/np.max(Y)


def hcosh_edge(tlist: np.ndarray,
               t_delay: float,
               t_width: float,
               t_plateau: float):
    
    zero1, zero2 = np.log(2-np.sqrt(3)), np.log(2+np.sqrt(3))
    t_width = max(1e-1, t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    w = t_width
    d = t_delay
    p = t_plateau

    x1 = zero1 * ((tlist - d)/(0.5 * w) - 1)
    x2= zero2 * ((tlist - d - w * 0.5 - p)/(0.5 * w))

    Y += ( tlist < d ) * 0
    Y += ( (tlist >= d) & (tlist <= d + w/2) ) * ( - np.cosh(x1) + 2)
    Y += ((tlist > d + w/2) & (tlist < d + w/2 + p)) * 1
    Y += ((tlist >= d + w/2 + p) & (tlist <= d + w + p)) * ( - np.cosh(x2) + 2)
    Y += ( tlist > (d + w + p) ) * 0

    return Y


def tanh_edge(tlist: np.ndarray,
              t_delay: float,
              t_width: float,
              t_plateau: float,
              epsilon: float):
    
    t_width = max(1e-1, t_width)
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


def hyper_edge(tlist: np.ndarray,
               t_delay: float,
               t_width: float,
               t_plateau: float,
               epsilon: float):
    
    t_width = max(1e-1, t_width)
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


def linear_edge(tlist: np.ndarray,
                t_delay: float,
                t_width: float,
                t_plateau: float):
    
    t_width = max(1e-1, t_width)
    Y = np.zeros(len(tlist) , dtype=np.complex128)
    w = t_width
    d = t_delay
    p = t_plateau
    tau = w/2

    ramp_slope = 1/tau

    Y += ( tlist < d ) * 0
    Y += ( (tlist >= d) & (tlist <= d + w/2) ) * ramp_slope * (tlist - d)
    Y += ((tlist > d + w/2) & (tlist < d + w/2 + p)) * 1
    Y += ( (tlist >= d + w/2 + p) & (tlist <= d + w + p) ) * -ramp_slope * (tlist - (d + w + p))
    Y += ( tlist > (d + w + p) ) * 0

    return Y


def linear_ramp_edge(tlist: np.ndarray,
                    t_delay: float,
                    t_width: float,
                    t_plateau: float,
                    cntrl: Literal['l', 'r'] = 'l'):
    
    t_width = max(1e-1, t_width)
    w = t_width
    d = t_delay
    p = t_plateau
    dt = tlist[1]-tlist[0]

    if cntrl == 'l':
        Y = linear_edge(tlist, d, w, p)
        Y *= 1 - ( (tlist >= d + w/2 + p) & (tlist <= d + w + p) ) 
    elif cntrl == 'r':
        d -= w/2 + dt
        Y = linear_edge(tlist, d, w, p)
        Y *= 1 - ( (tlist >= d) & (tlist <= d + w/2) ) 

    return Y


def cos_ramp_edge(tlist: np.ndarray,
                t_delay: float,
                t_width: float,
                t_plateau: float,
                cntrl: Literal['l', 'r'] = 'l'):
    
    t_width = max(1e-1, t_width)
    w = t_width
    d = t_delay
    p = t_plateau

    if cntrl == 'l':
        Y = cosine_edge(tlist, d, w, p)
        Y *= 1 - ( (tlist >= d + w/2 + p) & (tlist <= d + w + p) ) 
    elif cntrl == 'r':
        d -= w/2
        Y = cosine_edge(tlist, d, w, p)
        Y *= 1 - ( (tlist >= d) & (tlist <= d + w/2) ) 

    return Y

# #%%

# import numpy as np
# import matplotlib.pyplot as plt
# from typing import Literal


# def square_edge(tlist: np.ndarray,
#                 t_delay: float,
#                 t_plateau: float):
#     d = t_delay
#     p = t_plateau

#     return ( ((tlist >= d) & (tlist <= d + p) ))


# tlist = np.linspace(0, 100, 1000)
# t_delay = 10
# t_width = 40
# t_plateau = 20
# plt.plot(tlist, square_edge(tlist, t_delay, t_plateau))
