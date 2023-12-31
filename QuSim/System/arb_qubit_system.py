# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
from qutip import *
import QuSim.PulseGen.pulse_waveform as pw
import QuSim.Instruments.tools as tool

def hermitian_conjugate(matrix):
    # Convert the input list to a NumPy array
    np_matrix = np.array(matrix, dtype=np.complex_)
    
    # Compute the Hermitian conjugate
    hermitian_matrix = np_matrix.conj().T
    
    # Convert the result back to a list
    result = hermitian_matrix.tolist()
    
    return result
    
def cal_tensor(mylist):
        return tensor(*mylist)

class arb_qubit_system:
    """
    A class of a multiple qubit interacting system
    """

    def __init__(self, freq_list, inter_list, r=0, gamma_list=None, driving_list = None, bias_list = None):
        self.freq_list = freq_list # Energy of each state, each qubit
        self.num_q = len(self.freq_list) # Num of qubit
        self.inter_list = inter_list 
        if self.num_q > 1: self.r = r
        else: self.r = 0
        self.g = self.r # Coupling strength
        self.gamma_list = gamma_list # Lindblad operator

        self.q_dim_list = [len(freq) for freq in self.freq_list]
        self.total_dim = np.prod(np.array(self.q_dim_list))
        self.qeye_list = [qeye(n) for n in self.q_dim_list]
        self.H_0 = self.get_H_0()
        self.H_inter = self.get_H_inter()
        self.H = 2 * np.pi * (self.H_0 + self.H_inter)
        
        self.driving_list = driving_list
        self.bias_list = bias_list

        self.N = self.num_q * max(self.q_dim_list) + 1
        self.state_dic = enr_state_dictionaries(self.q_dim_list, self.N)

    def get_H_0(self):
        H_0 = 0
        for index, freq in enumerate(self.freq_list):
            qeye_list = self.qeye_list.copy()
            dummy_list = [] 
            empty_list = [0 for i in range(self.q_dim_list[index])]
            for i, w in enumerate(freq):
                row = list(np.copy(empty_list))
                row[i] = w
                dummy_list.append(row)
            dummy_qobj = Qobj(dummy_list)
            qeye_list[index] = dummy_qobj
            H_0 += cal_tensor(qeye_list)
        H_0.dims = [[self.total_dim], [self.total_dim]]
        return H_0

    def get_H_inter(self):
        H_inter = 0
        if self.g == 0: return 0
        if self.num_q > 1:
            for q_index1 in range(self.num_q - 1):
                for q_index2 in range(q_index1 + 1, self.num_q):
                    qeye_list = self.qeye_list.copy()
                    qeye_list[q_index1] = self.get_V(q_index1)
                    qeye_list[q_index2] = self.get_V(q_index2)
                    H_inter += self.g[q_index1][q_index2] * cal_tensor(qeye_list)
        H_inter.dims = [[self.total_dim], [self.total_dim]]
        return H_inter

    def get_V(self, qubit_index):
        q_dim = self.q_dim_list[qubit_index]
        empty_list = [0 for i in range(q_dim)]
        coupling_dic = self.inter_list[qubit_index]
        V_upper = []
        for i in range(q_dim - 1):
            row = list(np.copy(empty_list))
            for j in range(i+1, q_dim):
                row[j] = tool.get_v_element(coupling_dic, i, j)
            V_upper.append(row)
        # for i, v in enumerate(coupling):
        #     row = list(np.copy(empty_list))
        #     row[i+1] = v
        #     V_upper.append(row)
        V_upper.append(empty_list)
        # Get hermitian conjugate
        V_lower = hermitian_conjugate(V_upper)
        V_matrix = Qobj(list(np.array(V_upper) + np.array(V_lower)))
        # print(V_matrix)
        return V_matrix
    
    def get_H_XY_drive(self, qubit_index):
        if self.driving_list is None: return 0

        q_dim = self.q_dim_list[qubit_index]
        driving_dic = self.driving_list[qubit_index]
        empty_list = [0 for i in range(q_dim)]
        H_XY_upper = []
        for i in range(q_dim -1 ):
            row = list(np.copy(empty_list))
            for j in range(i+1, q_dim):
                row[j] = tool.get_XY_element(driving_dic, i, j)
            H_XY_upper.append(row)
        H_XY_upper.append(empty_list)
        H_XY_lower = hermitian_conjugate(H_XY_upper)
        H_XY_matrix = Qobj(list(np.array(H_XY_upper) + np.array(H_XY_lower)))
        return H_XY_matrix
    
    def get_H_Z_bias(self, qubit_index):
        q_dim = self.q_dim_list[qubit_index]
        a = destroy(q_dim)
        a_dagger = dag(a)
        if self.bias_list is None:
            return a_dagger * a
        # Add interface here...
        empty_list = [0 for i in range(q_dim)]
        H_Z = []
        H_Z.append(empty_list)
        for i in range(q_dim - 1):
            row = list(np.copy(empty_list))
            row[i+1] = tool.get_Z_element(self.bias_list[qubit_index], i+1)
            H_Z.append(row)
        H_Z_matrix = Qobj(H_Z)
        return H_Z_matrix

    def get_state_index(self, n, freq_threshold = 1e-6, deg_threshold = 5e-3, deg_round = 7):
        '''
        n: tuple
            e.g., (0,0,0), (1,0,1)
        
        freq_threshold: float, double
            The threshold of qubit frequency difference that 
            is recognized as degeneracies happening
        
        deg_threshold: float, double
            The threshold of probability amplitude that between
            superposition states constructed by energy
            degenerated states
        
        deg_round: int
            The number of decimal number that will be rounded in 
            estimating the probability amplitude of each degenerated 
            state.
        '''
        state_index = self.state_dic[1][n]
        eigen_list = [np.abs(arr[state_index][0][0]) for arr in self.H.eigenstates()[1]]

        max_value = max(eigen_list)
        max_index = np.argmax(np.array(eigen_list))
        # print(len(eigen_list))
        # print('eigen_list = {}'.format(eigen_list))
        # print('max_value = {}',format(max_value))
        # print('max_index = {}'.format(max_index))

        sim_index = tool.find_similar_indices(np.array(self.w), freq_threshold)
        if len(sim_index) > 0: 
            degenerate_index = []
            for index, value in enumerate(eigen_list):
                if index == max_index: continue
                if np.abs(value - max_value) < deg_threshold: # Not sufficient to say degeneracy is appears
                    prob_amp_list = []

                    # print(value)
                    # Exam the state 
                    for row in self.H.eigenstates()[1][index]:
                        prob_amp_list.append(np.round(np.abs(row[0][0]), deg_round))

                    # Count the number of equal array elements
                    deg_prob_amp_list, num_degen_list = np.unique(np.array(prob_amp_list), return_counts=True)
                    
                    # print(deg_prob_amp_list)
                    # print(num_degen_list)
                    # Extracting the maximum degenerate
                    i_max = np.argmax(deg_prob_amp_list)
                    num_degen = num_degen_list[i_max]
                    deg_prob_amp = deg_prob_amp_list[i_max]
                    
                    # print(deg_prob_amp)
                    # print(num_degen)

                    if num_degen > 1:
                        degenerate_index.append(index)

            # print(degenerate_index)
            if len(degenerate_index) > 0:
                degenerate_index.append(max_index)
                deg_index_arr = np.sort(np.array(degenerate_index))

                # Effective excitation number
                num_excit = 0
                for wi in sim_index:
                    num_excit += n[wi]

                count_n = 0
                for ii,wi in enumerate(sim_index):
                    if n[wi] != 0:
                        count_n += 1 * 2**(ii)
                max_index = deg_index_arr[count_n-1]
        # print(max_index)

        return max_index
    
    def get_eigenstates_energy(self, n, freq_threshold = 1e-6, deg_threshold = 5e-3, deg_round = 7):
        '''
        n: tuple
            e.g., (0,0,0), (1,0,1)
        
        freq_threshold: float, double
            The threshold of qubit frequency difference that 
            is recognized as degeneracies happening
        
        deg_threshold: float, double
            The threshold of probability amplitude that between
            superposition states constructed by energy
            degenerated states
        
        deg_round: int
            The number of decimal number that will be rounded in 
            estimating the probability amplitude of each degenerated 
            state.
        '''
        max_index = self.get_state_index(n, freq_threshold = 1e-6, deg_threshold = 5e-3, deg_round = 7)
        eigen_val_state = self.H.eigenstates()
        eigenstates, eigenenergies = eigen_val_state[1][max_index], eigen_val_state[0][max_index]/2/np.pi
        
        # Return a qobj eigenstate, energy level magnitude, and the index of the energy  level
        return eigenstates, eigenenergies.real, max_index
    
    def co_list(self):
        co_list = []
        if self.gamma_list is None:
            return co_list
        for q_index in range(self.num_q):
            qeye_list = self.qeye_list.copy()
            # Define creation/annihilation operator to define the collapse operator
            a = destroy(self.q_dim_list[q_index])
            a_dagger = dag(a)
            gamma_sum = 0
            if "up" in self.gamma_list[q_index]: 
                if self.gamma_list[q_index]["up"] != 0:
                    gamma_sum += np.sqrt(self.gamma_list[q_index]["up"]) * a_dagger
            if "down" in self.gamma_list[q_index]:
                if self.gamma_list[q_index]["down"] != 0:
                    gamma_sum += np.sqrt(self.gamma_list[q_index]["down"]) * a
            if "z" in self.gamma_list[q_index]:
                if self.gamma_list[q_index]["z"] != 0:
                    gamma_sum += np.sqrt(self.gamma_list[q_index]["z"] / 2) * a_dagger * a
            # Append the single qubit collapse operator to the collapse operator list
            if gamma_sum != 0:
                qeye_list[q_index] = gamma_sum
                collapse_operator = cal_tensor(qeye_list)
                collapse_operator.dims = [[self.total_dim], [self.total_dim]]
                co_list.append(collapse_operator)
        return co_list
    
    def system_dynamics_mesolve(self, simulation_option, pulse_sequence, option = Options(rtol=1e-8)):
        """
        A method
        
        """
        state_list = simulation_option["initial_state"]
        result_list, angle_list = [], []
        for state in state_list:
            H_d = []
            H_d.append(self.H)
            for pulse in pulse_sequence:
                H_d.append(self.send_pulse(pulse, simulation_option))
            initial_state = state
            simulation_step = simulation_option["simulation_step"]
            simulation_time = simulation_option["simulation_time"]
            # Set up master equation solver
            result, angle = self.master_eq_solver(H_d, simulation_time, simulation_step, initial_state, option = Options(rtol=1e-8))
            result_list.append(result)
            angle_list.append(angle)
        return result_list, angle_list
    
    # def system_dynamics_propagator(self, simulation_option, pulse_sequence, option = Options(rtol=1e-8)):
    #     H_d = []
    #     H_d.append(self.H)
    #     for pulse in pulse_sequence:
    #         H_d.append(self.send_pulse(pulse, simulation_option))
    #     simulation_step = simulation_option["simulation_step"]
    #     simulation_time = simulation_option["simulation_time"]
    #     tlist=np.linspace(0, simulation_time, simulation_step)
    #     result = propagator(H_d, tlist, self.co_list(), {} , option)
    #     return result
    
    def system_dynamics_propagator(self, simulation_option, pulse_sequence, option = Options(rtol=1e-8)):
        H_d = []
        H_d.append(self.H)
        for pulse in pulse_sequence:
            H_d.append(self.send_pulse(pulse, simulation_option))
        simulation_step = simulation_option["simulation_step"]
        simulation_time = simulation_option["simulation_time"]
        tlist=np.linspace(0, simulation_time, simulation_step)
        # tlist = simulation_time
        result = propagator(H_d, tlist, self.co_list(), {} , option, parallel=True, progress_bar=True)
        return result
    
    def send_pulse(self, pulse, simulation_option):
        """
        Construct dynamical component of the Hamiltonian H_d
        """
        if pulse['q_index'] > self.num_q - 1: ValueError('Invalid qubit index:'+ pulse['pulse_index']+ ', q_index = ' + pulse['q_index'])
        pulse["amplitude"] *= 2 * np.pi
        if pulse["type"] == "XY":
            H_drive = self.H_XY_drive(pulse, simulation_option)
        elif pulse["type"] == "Z":
            H_drive = self.H_Z_bias(pulse, simulation_option)
        else:
            raise ValueError("Invalid pulse type:"+ pulse['pulse_index']+ ', q_index = ' + pulse['pulse_index']+ ', type = '  + pulse["type"])
        pulse["amplitude"] /= 2 * np.pi
        return H_drive

    def H_XY_drive(self, pulse, simulation_option):
        """
        Define the Hamiltonian of the system under XY pulse driving

        t_width: float
            The sum of widths of the rising and lowering edges 
        t_plateau: float
            The width of plateau
        t_delay:
            The delay of the starting point of the pulse 
        pulse_shape:
            The waveform of the envelope
        ampli: float
            Amplitude of the XY pulse envelope
        freq: float
            Frequency of the XY pulse carrier
        phase: float
            Phase of the carrier
        q_index:
            The index of qubit that applied pulse to

        ===========================

        More features:
            - More pulse shape
            T.B.C.
        """
        q_index = pulse["q_index"]
        # Get pulse
        pulse_lib_class = pw.pulse_lib(pulse)
        XY_pulse = pulse_lib_class.get_pulse(simulation_option)
        # Get self driving Hamiltonian
        H_XY_drive = self.get_H_XY_drive(q_index)
        # Tensor product, get the total driving Hamiltonian
        qeye_list = self.qeye_list.copy()
        qeye_list[q_index] = H_XY_drive
        H_XY = cal_tensor(qeye_list)
        H_XY.dims = [[self.total_dim], [self.total_dim]]
        return [H_XY, XY_pulse]

    def H_Z_bias(self, pulse, simulation_option):
        """
        Define the Hamiltonian of the system under Z pulse biasing

        t_width: float
            The sum of widths of the rising and lowering edges 
        t_plateau: float
            The width of plateau
        t_delay:
            The delay of the starting point of the pulse 
        pulse_shape:
            The waveform of the envelope
        ampli: float
            Amplitude of the Z pulse
        q_index:
            The index of qubit that applied pulse to
        """
        q_index = pulse["q_index"]
        # Get flux pulse
        pulse_lib_class = pw.pulse_lib(pulse)
        flux_pulse = pulse_lib_class.get_pulse(simulation_option)

        # Get self bias Hamiltonian
        H_Z_bias = self.get_H_Z_bias(q_index)
        # Tensor product, get the total bias Hamiltonian
        qeye_list = self.qeye_list.copy()
        qeye_list[q_index] = H_Z_bias
        H_Z = cal_tensor(qeye_list)
        H_Z.dims = [[self.total_dim], [self.total_dim]]
        return [H_Z, flux_pulse]

    def master_eq_solver(self, H_d, t_simulation, simulation_step, initial_state, option = Options(rtol=1e-8)):
        tlist=np.linspace(0, t_simulation, simulation_step)
        result = mesolve(H_d, initial_state, tlist, self.co_list(), [], options=option) 
        angle = np.angle(initial_state.dag() * result.states[-1])
        return result, angle

    def get_data_list(self, result_list, simulation_option, state_list):
        simulation_step = simulation_option["simulation_step"]
        data_list = []
        index = 0
        for state in state_list:
            data_list_dummy = []
            for ii in range(0, simulation_step):
                data_list_dummy.append(np.abs(((result_list.states[ii]).dag()*state)[0][0][0])**2)
            data_list.append(data_list_dummy)
            index += 1
        return data_list
    
    def get_data_list_density(self, result, simulation_option, state_list):
        simulation_step = simulation_option["simulation_step"]
        data_list = []
        index = 0
        for state in state_list:
            data_list_dummy = []
            for ii in range(0, simulation_step):
                data_list_dummy.append(((result.states[ii] * result[index].states[ii].dag()) * (state * state.dag())).tr())
            data_list.append(data_list_dummy)
            index += 1
        return data_list
