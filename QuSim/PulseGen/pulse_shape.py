# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import qusim.PulseGen.edges as edges
import numpy as np

def cosine(tlist, pulse):
    # If the width of the pulse is small, then return a square pulse.
    if np.abs(pulse['t_width']) < 1e-5: return square(tlist, pulse)
    # print('pulse_shape.py, amp = {}'.format(pulse['amplitude']))
    return pulse['amplitude'] * edges.cosine_edge(tlist, pulse['t_delay'], pulse['t_width'], pulse['t_plateau'])

def hcosine(tlist, pulse): # Half cosine
    if np.abs(pulse['t_width']) < 1e-5: return square(tlist, pulse)
    a = 0.3; b = 1;
    return pulse['amplitude'] * edges.cosine_edge(tlist, pulse['t_delay'], pulse['t_width'], pulse['t_plateau'], a, b)

def cosh(tlist, pulse):
    if np.abs(pulse['t_width']) < 1e-5: return square(tlist, pulse)
    return pulse['amplitude'] * edges.hcosh_edge(tlist, pulse['t_delay'], pulse['t_width'], pulse['t_plateau'])

def square(tlist, pulse):
    return pulse['amplitude'] * ( ((tlist >= pulse['t_delay']) & (tlist <= pulse['t_delay'] + pulse['t_width'] + pulse['t_plateau']) ) + 0j )

def tanh(tlist, pulse):
    if 'epsilon' in pulse:
        return pulse['amplitude'] * edges.tanh_edge(tlist, pulse['t_delay'], pulse['t_width'], pulse['t_plateau'], pulse['epsilon'])
    else: raise ValueError("Missing variable epsilon: pulse_index = " + str(pulse['pulse_index']))

def hyper(tlist, pulse):
    if 'epsilon' in pulse:
        return pulse['amplitude'] * edges.hyper_edge(tlist, pulse['t_delay'], pulse['t_width'], pulse['t_plateau'], pulse['epsilon'])
    else: raise ValueError("Missing variable epsilon: pulse_index = " + str(pulse['pulse_index']))

PULSE_SHAPE_SET = {
    'square': square, # Square pulse
    'cosine': cosine, # Cosine
    'hcosine': hcosine, # Half cosine
    'cosh': cosh, # Cosh function with parts greater than zero.
    'tanh': tanh,
    'hyper': hyper # Hyperbolic cosh share, with control parameter \epsilon
}

PULSE_SHAPE = PULSE_SHAPE_SET.keys()