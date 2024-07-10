# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
from scipy.signal import find_peaks

def find_union(*lists):
    # Convert each list to a set
    sets = [set(lst) for lst in lists]
    # Find the intersection of all sets
    intersection_set = set.union(*sets)
    # Convert the intersection set back to a list
    intersection_list = list(intersection_set)
    return intersection_list

class SwapAvoidCrossing:
    
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
    
