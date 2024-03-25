# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import os
import pickle
import numpy as np
from scipy.signal import find_peaks

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

def get_XY_element(driving_dic, index1, index2):
    if index1 >= index2: raise ValueError("Invalid index: index1 >= index2 when getting the H_XY_drive")
    key = f"W{index1}{index2}"
    if key in driving_dic: return driving_dic[key]
    else: return 0

def get_Z_element(bias_dic, index):
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

# class swap_avoid_crossing_old:
#     def __init__(self, energy_level_list, w_scan, threshold):
#         self.energy_level = np.array(energy_level_list)
#         self.energy_level_T = self.energy_level.T
#         self.w_scan = w_scan
#         self.threshold = threshold
#         self.last_i = 50

#     def swap(self):
#         num_level = len(self.energy_level_T)
#         num_step = len(self.w_scan)
#         last_indices = [[] for i in range(self.last_i)]
#         for i in range(num_step):
#             E_diff = []
#             for k in range(num_level-1):
#                 E_diff.append(self.energy_level_T.T[i][k+1] - self.energy_level_T.T[i][k])
#             indices = self.find_indices_greater_than_threshold(E_diff)
#             indices, last_indices = self.exam_indices(indices, last_indices)
#             if indices == []: pass;
#             else:
#                 print("yes, swap, state indices = {}, w = {}".format(indices, self.w_scan[i]))
#                 # print("last_two_i = {}".format(last_indices))
#                 for index in indices:
#                     self.energy_level_T = self.swap_elements(self.energy_level_T, index, i)
#         return self.energy_level_T.T
    
#     def swap_elements(self, my_list, index1, index2):
#         # Check if the given index is valid
#         if index1 >= len(my_list) - 1 or index2 > len(my_list[0]) - 1:
#             raise ValueError("Invalid index")
#         my_list_copy = np.copy(my_list)
#         # Swap the elements after index i in the two sublists
#         for j in range(index2, len(my_list[0])):
#             my_list[index1][j], my_list[index1 + 1][j] = my_list_copy[index1 + 1][j], my_list_copy[index1][j]
#         return np.array(my_list)
    
#     def find_indices_greater_than_threshold(self, E_diff):
#         indices = [i for i, val in enumerate(E_diff) if np.abs(val) < self.threshold]
#         return indices
    
#     def exam_indices(self, current_indices, last_indices):
#         last_indices_cp = list(np.copy(last_indices))
#         last_indices = find_union(*last_indices)
#         pop_list = []
#         # print("========================")
#         # print("last_indices union = {}".format(last_indices))
#         # print("original indices,{}".format(current_indices))
#         for k, index in enumerate(np.copy(current_indices)):
#             if index in last_indices:
#                 pop_list.append(k)
#         for k in sorted(pop_list, reverse=True):
#             current_indices.pop(k)
#         # print("pop current indices,{}".format(current_indices))
#         last_indices_cp.insert(0, current_indices)
#         last_indices_cp.pop()
#         # print("========================")
#         return current_indices, last_indices_cp
        
class swap_avoid_crossing:
    
    def __init__(self, energy_level, w_scan, threshold):
        self.energy_level = np.array(energy_level)
        self.energy_level_T = self.energy_level.T
        self.w_scan = w_scan
        self.threshold = threshold
        self.last_i = 50 #
        self.E_diff = self.get_E_diff()
        self.minima_points = self.find_avoid_crossing_point()
        self.valid_swap_points = self.valid_swap_exam()
        self.energy_level_T_swap = np.copy(self.energy_level_T)
    
    def get_E_diff(self):
        num_level = len(self.energy_level_T)
        num_step = len(self.w_scan)
        E_diff = np.empty([num_level-1, num_step])
        for i in range(num_level-1):
            for j in range(num_step):
                E_diff[i][j] = np.abs(self.energy_level_T[i+1][j] - self.energy_level_T[i][j])
        return E_diff

    def find_avoid_crossing_point(self):
        minima_points = []
        for i,differ in enumerate(self.E_diff):
            minima_indices = find_peaks(-np.array(differ))[0]
            minima_values = np.array(differ)[minima_indices]
            for j,case in enumerate(minima_indices):
                combined_infor = [case, minima_values[j], i, i+1]
                minima_points.append(combined_infor)
        # Sort the minima_points list base on the minima_indices from low to high
        sorted_minima_points = sorted(minima_points, key=lambda x: x[0])
        return sorted_minima_points

    def valid_swap_exam(self):
        valid_swap_points = []
        for i,point in enumerate(self.minima_points):
            if np.abs(point[1]) < self.threshold:
                valid_swap_points.append(point)
        for i,point in enumerate(valid_swap_points):
            valid_swap_points = self.update_point(i, point[2], point[3],valid_swap_points)
        return valid_swap_points
    
    def update_point(self, index, index_low, index_high, valid_swap_points):
        # print(valid_swap_points[index+1:])
        for i,point in enumerate(valid_swap_points[index+1:]):
            if point[2] == index_low:
                valid_swap_points[index+1 + i][2] = index_high
            elif point[2] == index_high:
                valid_swap_points[index+1 + i][2] = index_low
            if point[3] == index_low:
                valid_swap_points[index+1 + i][3] = index_high
            elif point[3] == index_high:
                valid_swap_points[index+1 + i][3] = index_low
        return valid_swap_points

    def swap(self):
        for point in self.valid_swap_points:
            # print(point)
            self.energy_level_T_swap = self.swap_elements(self.energy_level_T_swap, point[2], point[3], point[0])
        return self.energy_level_T_swap.T
    
    def swap_elements(self, my_list, index1_low, index1_high, index2):
        # Check if the given index is valid
        if index1_low >= len(my_list) - 1 or index2 > len(my_list[0]) - 1:
            raise ValueError("Invalid index")
        my_list_copy = np.copy(my_list)
        # Swap the elements after index i in the two sublists
        for j in range(index2, len(my_list[0])):
            my_list[index1_low][j], my_list[index1_high][j] = my_list_copy[index1_high][j], my_list_copy[index1_low][j]
        return np.array(my_list)

def find_similar_indices(arr, threshold):
    similar_indices = []
    
    for i in range(len(arr)):
        for j in range(i+1, len(arr)):
            if np.abs(arr[j] - arr[i]) < threshold:
                similar_indices.extend([i, j])

    return list(set(similar_indices))