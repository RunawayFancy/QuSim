import numpy as np
from qutip import *
import pulse_waveform as pw

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
    
    g: array like
        Two dimensional array, encoding
        the coupling strength between each
        two qubit
    """

    def __init__(self, N, q_dim, w, alpha, r, gamma_list):
        self.N = N
        self.w = w
        self.num_q = len(self.w)  # Number of qubits
        self.alpha = alpha
        if self.num_q > 1: self.r = r
        else: self.r = 0
        self.g = self.r
        self.q_dim = q_dim
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

    def system_dynamics_mesolve(self, simulation_option, pulse_sequence):
        state_list = simulation_option["initial_state"]
        result_list, angle_list = [], []
        for state in state_list:
            H_d = []
            H_d.append(self.H)
            for pulse in pulse_sequence:
                H_d.append(self.send_pulse(pulse))
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
            H_d.append(self.send_pulse(pulse))
        simulation_step = simulation_option["simulation_step"]
        simulation_time = simulation_option["simulation_time"]
        tlist=np.linspace(0, simulation_time, simulation_step)
        option = Options(rtol=1e-8)
        result = propagator(H_d, tlist, self.co_list(), {} , option)
        return result

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
        # if self.num_q == 2: H_inter += self.g[0][0]
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
    
    # def get_eigenstates_energy(self, n):
    #     '''
    #     n: tuple
    #         e.g., (0,0,0), (1,0,1)
    #     '''
    #     max_index = self.state_dic[1][n]
    #     eigen_val_state = self.H.eigenstates()
    #     eigenstates, eigenenergies = eigen_val_state[1][max_index], eigen_val_state[0][max_index]/2/np.pi
    #     # Return a qobj eigenstate, energy level magnitude, and the index of the energy  level
    #     return eigenstates, eigenenergies, max_index
    
    def send_pulse(self, pulse):
        if pulse['q_index'] > self.num_q - 1: ValueError('Invalid qubit index:'+ pulse['pulse_index']+ ', q_index = ' + pulse['q_index'])
        pulse["amplitude"] *= 2 * np.pi
        if pulse["type"] == "XY":
            H_drive = self.H_XY_drive(pulse)
        elif pulse["type"] == "Z":
            H_drive = self.H_Z_bias(pulse)
        else:
            raise ValueError("Invalid pulse type:"+ pulse['pulse_index']+ ', q_index = ' + pulse['pulse_index']+ ', type = '  + pulse["type"])
        pulse["amplitude"] /= 2 * np.pi
        return H_drive

    def H_XY_drive(self, pulse):
        """
        Define the Hamiltonian of the system under XY pulse driving
        t_g: float
            The gate time
        t_rising: float
            The rising and lowering time of XY pulse
        ampli: float
            Amplitude of the XY pulse
        freq: float
            Frequency of the XY pulse
        
        ===========================

        More features:
            - DRAG pulse regime
            - More pulse shape
            T.B.C.
        """
        q_index = pulse["q_index"]
        # Get pulse
        pulse_lib_class = pw.pulse_lib(pulse)
        XY_pulse = pulse_lib_class.get_pulse()
        
        return [self.a_dagger_list[q_index] * self.a_list[q_index], XY_pulse]

    def H_Z_bias(self, pulse):
        q_index = pulse["q_index"]
        # Get flux pulse
        pulse_lib_class = pw.pulse_lib(pulse)
        flux_pulse = pulse_lib_class.get_pulse()
        return [self.a_list[q_index] + self.a_dagger_list[q_index], flux_pulse]
    
    def co_list(self):
        co_list = []
        if self.gamma_list is None:
            return co_list
        for q_index in range(self.num_q):
            # Get collapse up operator
            co_list.append(np.sqrt(self.gamma_list[q_index]["up"]) * self.a_dagger_list[q_index])
            # Get collapse down operator
            co_list.append(np.sqrt(self.gamma_list[q_index]["down"]) * self.a_list[q_index])
            # Get collapse z operator
            # Question marks: L_z = sqrt(2 Gamma_Z) a^dagger a
            co_list.append(np.sqrt(2 * self.gamma_list[q_index]["z"]) * self.a_dagger_list[q_index] * self.a_list[q_index])

        return co_list

    def master_eq_solver(self, H_d, t_simulation, simulation_step, initial_state):
        tlist=np.linspace(0, t_simulation, simulation_step)
        option = Options(rtol=1e-8)
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