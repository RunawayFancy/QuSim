# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import os
import pickle
import numpy as np
import itertools

# SwapAvoidCrossing is defined in avoid_crossing.py; re-export it here so that
# `qusim.instruments.tools.SwapAvoidCrossing` keeps working after de-duplication.
from .avoid_crossing import SwapAvoidCrossing

def r2matrix(r_dic, frequency_list):
    n = len(frequency_list)
    r = [[0] * n for _ in range(n)]
    # if n ==1: r[0][0] = r_dic.get(f"r{1}{2}", 0)
    # else:
    for i in range(n):
        for j in range(n):
            if j >= i:
                key = f"r{i+1}{j+1}"
                if key in r_dic:
                    r[i][j] = r_dic.get(f"r{i+1}{j+1}", 0)
                else: r[i][j] = 0
    return r

def get_v_diag_element(v_dic, index1, index2):
    if index1 != index2: raise ValueError("Invalid index: index1 != index2 when getting the coupling strength")
    key = f"v{index1}{index2}"
    if key in v_dic: return v_dic[key]
    else: return 0

def get_v_element(v_dic, index1, index2):
    if index1 >= index2: raise ValueError("Invalid index: index1 >= index2 when getting the coupling strength")
    key = f"v{index1}{index2}"
    if key in v_dic:
        # print(key) 
        return v_dic[key]
    else: return 0

def get_xy_element(driving_dic, index1, index2):
    if index1 >= index2: raise ValueError("Invalid index: index1 >= index2 when getting the XY drive element")
    key = f"W{index1}{index2}"
    if key in driving_dic: return driving_dic[key]
    else: return 0

def get_z_element(bias_dic, index):
    key = f"Z{index}{index}"
    if key in bias_dic: return bias_dic[key]
    else: return 0

def write_data(output_file, data):
    # Write the simu_data dictionary to a local file using pickle
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "wb") as file:
        pickle.dump(data, file)

def load_data(input_file):
    # Load the data from the file using pickle
    with open(input_file, "rb") as file:
        loaded_data = pickle.load(file)
    return loaded_data

def grad(sh): 
  if(len(sh) <=1) :return 0*sh   ; 
  else : return np.gradient(sh) ; 


def find_similar_indices(arr, threshold):
    similar_indices = []
    
    for i in range(len(arr)):
        for j in range(i+1, len(arr)):
            if np.abs(arr[j] - arr[i]) < threshold:
                similar_indices.extend([i, j])

    return list(set(similar_indices))


def get_combinations(num):
    numbers = range(num)
    combinations = list(itertools.combinations(numbers, 2))
    
    return combinations