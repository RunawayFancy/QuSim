# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
from qusim.PulseGen.noise_config import * 
import numpy as np
from typing import Union

def noise_gen(
        noise_config: Union[GaussianNoiseConfig|RandomTeleNoiseConfig|JNNoiseConfig|OneOverFNoiseConfig],
        waveform: np.ndarray
    ) -> np.ndarray:
    """
    Generate noise to be added to a qubit control pulse.

    Returns:
    - noise: A numpy array containing the generated noise over time.
    """

    tlist = noise_config.noise_time_config.simopt.tlist
    if noise_config.noise_time_config.time_dependent:
        if noise_config.noise_time_config.tdbasefn:
            noise_base = noise_config.noise_time_config.tdbasefn(tlist, waveform)
        else:
            raise TypeError("Time dependent noise requires time-dependent base function `tdbasefn`.")
    else:
        noise_base = np.ones_like(tlist)

    if noise_config.type == '1/f':
        noise_arr = noise_base * noise_config.trigger(len(tlist))
    else:

        noise_segments = segmentize(tlist, noise_config.noise_time_config.tseg, noise_base)
        noise_multiplier = noise_config.trigger(len(noise_segments))

        noise_arr = np.array(np.concatenate([sublist * noise_multiplier[_i] for _i, sublist in enumerate(noise_segments)]))

    if noise_config.noise_time_config.tranfofn:
        noise_arr = noise_config.noise_time_config.tranfofn(noise_arr)

    noise_arr[(tlist < noise_config.noise_time_config.tstart) | (tlist > noise_config.noise_time_config.tstop)] = 0

    return noise_arr


def segmentize(tlist: np.ndarray, tseg: float, noise: np.ndarray) -> np.ndarray:
    segments = []
    start_idx = 0

    while start_idx < len(tlist):
        # Find the end index for the current segment
        end_idx = start_idx
        while end_idx < len(tlist) and (tlist[end_idx] - tlist[start_idx]) <= tseg:
            end_idx += 1
        
        # Add the segment to the list
        segments.append(noise[start_idx:end_idx])
        
        # Move to the next starting index
        start_idx = end_idx

    return np.array(segments, dtype='object')
