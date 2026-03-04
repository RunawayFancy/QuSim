# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
from qusim.PulseGen.noise_config import *
import numpy as np
from typing import Union, Optional, Tuple, Dict

# Cache for correlated random draws across channels
_CORR_RAND_CACHE: Dict[Tuple, np.ndarray] = {}

def _config_signature(noise_config) -> Tuple:
    ntype = noise_config.type
    if ntype == 'ga':
        return (ntype, noise_config.mean, noise_config.std, noise_config.amp)
    if ntype == 'rt':
        return (ntype, noise_config.high_val, noise_config.low_val, noise_config.switch_prob)
    if ntype == 'jn':
        return (ntype, noise_config.temperature, noise_config.resistance)
    if ntype == '1/f':
        return (ntype, noise_config.alpha, noise_config.scale, noise_config.cutoff_freq_low, noise_config.cutoff_freq_high)
    return (ntype,)

def _corr_key(noise_config, tlist: np.ndarray, rnd_len: int):
    corr_id = getattr(noise_config, "corr_id", None)
    if not corr_id:
        return None
    ntc = noise_config.noise_time_config
    return (
        corr_id,
        id(tlist),
        _config_signature(noise_config),
        rnd_len,
        ntc.tseg,
        ntc.tstart,
        ntc.tstop,
        ntc.time_dependent
    )

def _get_or_make_random(noise_config, tlist: np.ndarray, rnd_len: int) -> np.ndarray:
    key = _corr_key(noise_config, tlist, rnd_len)
    if key is not None and key in _CORR_RAND_CACHE:
        return _CORR_RAND_CACHE[key]
    rnd = noise_config.trigger(rnd_len)
    if key is not None:
        _CORR_RAND_CACHE[key] = rnd
    return rnd

def clear_corr_noise_cache(corr_id: Optional[str] = None) -> None:
    """
    Clear cached correlated random draws. If corr_id is None, clear all.
    """
    if corr_id is None:
        _CORR_RAND_CACHE.clear()
        return
    keys = [k for k in _CORR_RAND_CACHE.keys() if k and k[0] == corr_id]
    for k in keys:
        del _CORR_RAND_CACHE[k]

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
        rnd = _get_or_make_random(noise_config, tlist, len(tlist))
        noise_arr = noise_base * rnd
    else:
        noise_segments = segmentize(tlist, noise_config.noise_time_config.tseg, noise_base)
        rnd = _get_or_make_random(noise_config, tlist, len(noise_segments))
        noise_arr = np.array(
            np.concatenate([sublist * rnd[_i] for _i, sublist in enumerate(noise_segments)]),
            dtype='float64'
        )

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
