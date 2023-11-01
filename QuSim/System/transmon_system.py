# -*- coding: utf-8 -*-
"""
@author: Pan Shi, Jiheng Duan
"""
import numpy as np
from qutip import *
import QuSim.PulseGen.pulse_waveform as pw

class qubit_system:
    """
    A class of a multiple qubit (transmon) interacting system
    N: int
        Dimension of the Hilbert space of every qubit (number of energy level considered in the simulation)

    w: array like
        One dimensional array, encoding frequency of each qubit

    alpha: array like
        One dimensional array, encoding anharmonicity of each qubit
        For qubit, let alpha < 0
        For resonator, let alpha = 0
    
    r: array like
        Two dimensional array, encoding
        the coupling strength between each
        two qubit

    q_dim: array like
        One dimensional list, defined by
        `[qubit_dim for _ in range(len(w))]`
    
    gamma_list: array like
        One dimensional list, containing multiple dictionaries
        with numbers equal to the number of qubits.
    """

    def __init__(self, N, q_dim, w, alpha, r=0, gamma_list = None):
        self.w = w
        self.num_q = len(self.w)  # Number of qubits
        self.q_dim = q_dim # Qubit dimension
        if N != None: self.N = N # Turn on max photon num
        else: self.N = self.num_q * max(self.q_dim) + 1 # Turn off max photon num
        self.alpha = alpha
        if self.num_q > 1: self.r = r
        else: self.r = 0
        self.g = self.r
        self.a_list = self.get_a_list()  # Define the second quantization field operator
        self.a_dagger_list = [dag(a) for a in self.a_list]
        self.H_q, self.H_a = self.get_Hq_Ha()
        self.H_inter = self.get_H_inter()
        self.H = 2 * np.pi * (self.H_q + self.H_a + self.H_inter)
        self.gamma_list = gamma_list
        self.state_dic = enr_state_dictionaries(self.q_dim, self.N)
    
    def get_a_list(self):
        a_list = enr_destroy(self.q_dim, excitations=self.N)
        row_shape, col_shape = a_list[0].shape[0], a_list[0].shape[1]
        for a in a_list:
            a.dims = [[row_shape], [col_shape]]
        return a_list

    def get_Hq_Ha(self):
        """
        Calculate the qubit and anharmonicity Hamiltonian
        H_q: qubit Hamiltonian
        H_a: anharmonicity

        _q: variable for H_q
        _a: variable for anharmonicity
        """
        # Define the qubit and anharmo Hamiltonian
        H_q, H_a = 0, 0
        for q_index in range(self.num_q):
            H_q += self.w[q_index] * self.a_dagger_list[q_index] * self.a_list[q_index]
            H_a += self.alpha[q_index]/2 * self.a_dagger_list[q_index] * self.a_dagger_list[q_index] * self.a_list[q_index] * self.a_list[q_index]
        
        return H_q, H_a

    def get_H_inter(self):
        H_inter = 0
        if self.g == 0: return 0
        if self.num_q > 1:
            for q_index1 in range(self.num_q - 1):
                for q_index2 in range(q_index1 + 1, self.num_q): 
                    H_inter += self.g[q_index1][q_index2] * (self.a_list[q_index1] + self.a_dagger_list[q_index1]) * (self.a_list[q_index2] + self.a_dagger_list[q_index2])
        return H_inter

    def get_state_index(self, n):
        '''
        n: tuple
            e.g., (0,0,0), (1,0,1)
        '''
        state_index = self.state_dic[1][n]
        eigen_list = [np.abs(arr[state_index][0][0].real) for arr in self.H.eigenstates()[1]]

        max_value = max(eigen_list)
        for index, value in enumerate(eigen_list):
            if value == max_value:
                max_index = index
                break
        return max_index


    def get_eigenstates_energy(self, n):
        '''
        n: tuple
            e.g., (0,0,0), (1,0,1)
        '''
        max_index = self.get_state_index(n)
        eigen_val_state = self.H.eigenstates()
        eigenstates, eigenenergies = eigen_val_state[1][max_index], eigen_val_state[0][max_index]/2/np.pi
        
        # Return a qobj eigenstate, energy level magnitude, and the index of the energy  level
        return eigenstates, eigenenergies.real, max_index
    
    def co_list(self):
        co_list = []
        if self.gamma_list is None:
            return co_list
        for q_index in range(self.num_q):
            # Get collapse up operator
            if self.gamma_list[q_index]["up"] != 0:
                co_list.append(np.sqrt(self.gamma_list[q_index]["up"]) * self.a_dagger_list[q_index])
            # Get collapse down operator
            if self.gamma_list[q_index]["down"] != 0:
                co_list.append(np.sqrt(self.gamma_list[q_index]["down"]) * self.a_list[q_index])
            # Get collapse z operator
            # Question marks: L_z = sqrt(2 Gamma_Z) a^dagger a
            if self.gamma_list[q_index]["z"] != 0:
                co_list.append(np.sqrt(2 * self.gamma_list[q_index]["z"]) * self.a_dagger_list[q_index] * self.a_list[q_index])

        return co_list
    
    def system_dynamics_mesolve(self, simulation_option, pulse_sequence):
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
            result, angle = self.master_eq_solver(H_d, simulation_time, simulation_step, initial_state)
            result_list.append(result)
            angle_list.append(angle)
        return result_list, angle_list
    
    def system_dynamics_propagator(self, simulation_option, pulse_sequence):
        H_d = []
        H_d.append(self.H)
        for pulse in pulse_sequence:
            H_d.append(self.send_pulse(pulse, simulation_option))
        simulation_step = simulation_option["simulation_step"]
        simulation_time = simulation_option["simulation_time"]
        tlist=np.linspace(0, simulation_time, simulation_step)
        option = Options(rtol=1e-8)
        result = propagator(H_d, tlist, self.co_list(), {} , option)
        return result
    
    def send_pulse(self, pulse, simulation_option):
        """
        Construct dynamical component of the Hamiltonian H_d
        """
        if pulse['q_index'] > self.num_q - 1: ValueError('Invalid qubit index:'+ pulse['pulse_index']+ ', q_index = ' + pulse['q_index'])
        # print('send_pulse()_1, amp = {}'.format(pulse['amplitude']))
        pulse["amplitude"] *= 2 * np.pi
        if pulse["type"] == "XY":
            H_drive = self.H_XY_drive(pulse, simulation_option)
        elif pulse["type"] == "Z":
            H_drive = self.H_Z_bias(pulse, simulation_option)
        else:
            raise ValueError("Invalid pulse type:"+ pulse['pulse_index']+ ', q_index = ' + pulse['pulse_index']+ ', type = '  + pulse["type"])
        # print('send_pulse()_2, amp = {}'.format(pulse['amplitude']))
        pulse["amplitude"] /= 2 * np.pi
        # print('send_pulse()_3, amp = {}'.format(pulse['amplitude']))
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

        # print('H_XY_drive, amp = {}'.format(pulse['amplitude']))

        XY_pulse = pulse_lib_class.get_pulse(simulation_option)
        
        return [-1j*self.a_dagger_list[q_index] + 1j*self.a_list[q_index], XY_pulse]

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

        return [self.a_dagger_list[q_index] * self.a_list[q_index], flux_pulse]

    def master_eq_solver(self, H_d, t_simulation, simulation_step, initial_state):
        tlist = np.linspace(0, t_simulation, simulation_step)
        option = Options(rtol=1e-8)
        result = mesolve(H_d, initial_state, tlist, c_ops = self.co_list(), options = option) 
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