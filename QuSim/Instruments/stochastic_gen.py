# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
import copy

def gen_normal_rand_var(mean, std):
    """
    Generate a random variable from a Gaussian distribution.

    Parameters:
    - mean (float): Mean of the Gaussian distribution.
    - std (float): Standard deviation of the Gaussian distribution.

    Returns:
    - float: Random variable sampled from the Gaussian distribution.
    """
    return np.random.normal(mean, std)

def gen_poisson_rand_var(mean):
    """
    Generate a random variable from a Poisson distribution.

    Parameters:
    - mean (float): Mean of the Poisson distribution.

    Returns:
    - int: Random variable sampled from the Poisson distribution.
    """
    return np.random.poisson(mean)

def gen_rt_rand_var(p_switch, max_limit=1000):
    steps = 1
    while steps < max_limit:
        if np.random.rand() < p_switch:
            return steps
        else: 
            steps+=1

def find_closest_value_with_index(sorted_array, target):
    """
    Find the closest value to the target in a sorted array.

    Parameters:
    - sorted_array (list): A sorted list of numbers.
    - target (float): The target value.

    Returns:
    - tuple: A tuple containing the closest value and its index in the array.
    """
    low, high = 0, len(sorted_array) - 1

    while low <= high:
        mid = (low + high) // 2
        mid_value = sorted_array[mid]

        if mid_value == target:
            return mid_value, mid
        elif mid_value < target:
            low = mid + 1
        else:
            high = mid - 1

    # Check the values at low and high positions
    closest_low = sorted_array[low] if low < len(sorted_array) else None
    closest_high = sorted_array[high] if high >= 0 else None

    # Determine the closest value and its index
    if closest_low is None:
        return closest_high, high
    elif closest_high is None:
        return closest_low, low
    else:
        return (closest_low, low) if abs(target - closest_low) < abs(target - closest_high) else (closest_high, high)
    

def gen_stochastic_point(simulation_option, t_switching_mean):
    time_list = np.linspace(0, simulation_option["simulation_time"], simulation_option["simulation_step"])
    ti = np.copy(time_list[0])
    t = np.copy(ti)
    index = 0
    lenth = len(time_list)
    S = []
    index_list = []
    while t < time_list[-1]:
        interval = gen_poisson_rand_var(t_switching_mean)
        t += interval
        if t > time_list[-1]: break
        arr = np.array(time_list[index:lenth])
        t_c, index_new = find_closest_value_with_index(arr, t)
        S.append(t_c)
        index_list.append(index + index_new)
        index = np.copy(index_new)
    return S, index_list

def gen_rt_point(simulation_option, p_switch, max_limit = 1000):
    time_list = np.linspace(0, simulation_option["simulation_time"], simulation_option["simulation_step"])
    ti = np.copy(time_list[0])
    t = np.copy(ti)
    index = 0
    lenth = len(time_list)
    S = []
    index_list = []
    while t < time_list[-1]:
        interval = gen_rt_rand_var(p_switch, max_limit)
        t += interval
        if t > time_list[-1]: break
        arr = np.array(time_list[index:lenth])
        t_c, index_new = find_closest_value_with_index(arr, t)
        S.append(t_c)
        index_list.append(index + index_new)
        index = np.copy(index_new)
    return S, index_list

def gen_pulse_seq(S:list, ilist:list, dist_dic:dict, type:str, shape:str, q_indx:int):
    pulse_ex = {
        'pulse_index': 0,
        'type': type,
        'pulse_shape': shape,
        't_delay': 0, # unit in ns
        't_width': 0, # unit in ns
        't_plateau': 0, # unit in ns
        'freq': 0, # unit in GHz; Z pulse does not use it
        'phase': 0, # unit in rad; Z pulse does not use it
        'amplitude': 0, # XY: Rabi freq; Z: biased frequency
        'q_index': q_indx # 0, 1, 2 ...
    }
    pulse_sequence = []

    for ii, t_d in enumerate(S):
        pulse = copy.deepcopy(pulse_ex)

        pulse['pulse_index'] = ii + 1
        pulse['t_delay'] = t_d
        pulse['t_width'] = gen_normal_rand_var(dist_dic["t_width_mean"], dist_dic["t_width_std"])
        pulse['freq'] = gen_normal_rand_var(dist_dic["freq_mean"], dist_dic["freq_std"])
        pulse['phase'] = gen_normal_rand_var(dist_dic["phase_mean"], dist_dic["phase_std"])
        pulse['amplitude'] = gen_normal_rand_var(dist_dic["amplitude_mean"], dist_dic["amplitude_std"])

        pulse_sequence.append(pulse)

    return pulse_sequence
