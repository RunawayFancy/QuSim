import os
import pickle
import numpy as np

def r2matrix(r_dic):
    n = len(r_dic)
    r = [[0] * n for _ in range(n)]
    for i, key_i in enumerate(r_dic):
        for j, key_j in enumerate(r_dic):
            if j >= i:
                r[i][j] = r_dic.get(f"r{i+1}{j+1}", 0)
    return r

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

def find_union(*lists):
    # Convert each list to a set
    sets = [set(lst) for lst in lists]
    # Find the intersection of all sets
    intersection_set = set.union(*sets)
    # Convert the intersection set back to a list
    intersection_list = list(intersection_set)
    return intersection_list

class swap_avoid_crossing:
    def __init__(self, energy_level_list, w_scan, threshold):
        self.energy_level = np.array(energy_level_list)
        self.energy_level_T = self.energy_level.T
        self.w_scan = w_scan
        self.threshold = threshold
        self.last_i = 50

    def swap(self):
        num_level = len(self.energy_level_T)
        num_step = len(self.w_scan)
        last_indices = [[] for i in range(self.last_i)]
        for i in range(num_step):
            E_diff = []
            for k in range(num_level-1):
                E_diff.append(self.energy_level_T.T[i][k+1] - self.energy_level_T.T[i][k])
            indices = self.find_indices_greater_than_threshold(E_diff)
            indices, last_indices = self.exam_indices(indices, last_indices)
            if indices == []: pass;
            else:
                print("yes, swap, state indices = {}, w = {}".format(indices, self.w_scan[i]))
                # print("last_two_i = {}".format(last_indices))
                for index in indices:
                    self.energy_level_T = self.swap_elements(self.energy_level_T, index, i)
        return self.energy_level_T.T
    
    def swap_elements(self, my_list, index1, index2):
        # Check if the given index is valid
        if index1 >= len(my_list) - 1 or index2 > len(my_list[0]) - 1:
            raise ValueError("Invalid index")
        my_list_copy = np.copy(my_list)
        # Swap the elements after index i in the two sublists
        for j in range(index2, len(my_list[0])):
            my_list[index1][j], my_list[index1 + 1][j] = my_list_copy[index1 + 1][j], my_list_copy[index1][j]
        return np.array(my_list)
    
    def find_indices_greater_than_threshold(self, E_diff):
        indices = [i for i, val in enumerate(E_diff) if np.abs(val) < self.threshold]
        return indices
    
    def exam_indices(self, current_indices, last_indices):
        last_indices_cp = list(np.copy(last_indices))
        last_indices = find_union(*last_indices)
        pop_list = []
        # print("========================")
        # print("last_indices union = {}".format(last_indices))
        # print("original indices,{}".format(current_indices))
        for k, index in enumerate(np.copy(current_indices)):
            if index in last_indices:
                pop_list.append(k)
        for k in sorted(pop_list, reverse=True):
            current_indices.pop(k)
        # print("pop current indices,{}".format(current_indices))
        last_indices_cp.insert(0, current_indices)
        last_indices_cp.pop()
        # print("========================")
        return current_indices, last_indices_cp
        