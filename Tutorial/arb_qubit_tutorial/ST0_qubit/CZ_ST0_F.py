import sys
# This on is a path correcting code, which is used to find the path of QuSim.
sys.path.append('../../..')
import QuSim.DataPlot.plot_lib as pl
import QuSim.Instruments.tools as tools
import QuSim.Instruments.stochastic_gen as sg
import QuSim.System.arb_qubit_system as aqs
import QuSim.Instruments.angle as  ang
from QuSim.Instruments.angle import get_angle

import numpy as np
import matplotlib.pyplot as plt

from tqdm import *
from time import *
from sympy import*
import copy
import pickle
from qutip import*
import argparse

freq_list = [ # GHz
    [0, 10.35],
    [0, 10],
    [0, 10.25],
    [0, 10.15]
]

inter_list = [
    [{
        "v01": 0,
        "v00": 1,
        "v11": -1
    },
    {
        "v01": 0,
        "v00": 1,
        "v11": -1
    }
    ,
    {
        "v01": 0,
        "v00": 1,
        "v11": -1
    },
    {
        "v01": 0,
        "v00": 1,
        "v11": -1
    }],
    [{
        "v01": 1,
        "v00": 0,
        "v11": 0
    },
    {
        "v01": 1,
        "v00": 0,
        "v11": 0
    }
    ,
    {
        "v01": 1,
        "v00": 0,
        "v11": 0
    },
    {
        "v01": 1,
        "v00": 0,
        "v11": 0
    }]
    ,
    [{
        "v01": -1j,
        "v00": 0,
        "v11": 0
    },
    {
        "v01": -1j,
        "v00": 0,
        "v11": 0
    }
    ,
    {
        "v01": -1j,
        "v00": 0,
        "v11": 0
    },
    {
        "v01": -1j,
        "v00": 0,
        "v11": 0
    }]
]

r_dic = {
    "r12": 1e-3,
    "r23": 1e-3,
    "r34": 1e-3
}
r = tools.r2matrix(r_dic, freq_list) # Coupling strength

extra_list=None
## with pulse type XY
driving_list = [
    {
        "W01": 1j   # W01 represents \Omega_{01}, for the qubit with index 0
    },
    {
        "W01": 1j
    },
    {
        "W01": 1j
    },
    {
        "W01": 1j
    }
]

## with pulse type Z
bias_list = [ # Default, Z00 = 0
    {
        "Z11": 1
    },
    {
        "Z11": 1
    },
    {
        "Z11": 1
    },
    {
        "Z11": 1
    }
]

gamma_list = None

_system_arb1 = aqs.arb_qubit_system(freq_list, inter_list, r, extra_list, gamma_list, driving_list, bias_list)

state_0000, E_0000, _= _system_arb1.get_eigenstates_energy((0,0,0,0))

state_0001, E_0001, _= _system_arb1.get_eigenstates_energy((0,0,0,1))
state_0010, E_0010, _= _system_arb1.get_eigenstates_energy((0,0,1,0))
state_0100, E_0100, _= _system_arb1.get_eigenstates_energy((0,1,0,0))
state_1000, E_1000, _= _system_arb1.get_eigenstates_energy((1,0,0,0))

"""subspace"""
state_0101, E_0101, index_0101= _system_arb1.get_eigenstates_energy((0,1,0,1))
state_1001, E_1001, index_1001= _system_arb1.get_eigenstates_energy((1,0,0,1))
state_0110, E_0110, index_0110= _system_arb1.get_eigenstates_energy((0,1,1,0))
state_1010, E_1010, index_1010= _system_arb1.get_eigenstates_energy((1,0,1,0))

"""leakage"""
state_0011, E_0011, index_0011= _system_arb1.get_eigenstates_energy((0,0,1,1))
state_1100, E_1100, index_1100= _system_arb1.get_eigenstates_energy((1,1,0,0))

state_0111, E_0111, _= _system_arb1.get_eigenstates_energy((0,1,1,1))
state_1011, E_1011, _= _system_arb1.get_eigenstates_energy((1,0,1,1))
state_1101, E_1101, _= _system_arb1.get_eigenstates_energy((1,1,0,1))
state_1110, E_1110, _= _system_arb1.get_eigenstates_energy((1,1,1,0))

state_1111, E_1111, _= _system_arb1.get_eigenstates_energy((1,1,1,1))

state_ud = Qobj(np.array([0,0,1,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)
state_du = Qobj(np.array([0,1,0,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)

state_dudu = tensor(state_du,state_du)
state_duud = tensor(state_du,state_ud)
state_uddu = tensor(state_ud,state_du)
state_udud = tensor(state_ud,state_ud)

bisa_amp = 5.152e-3# 11.128e-3;
t_width = 35

cp_matix = pickle.load(open('cp_matrix.pkl', 'rb')) 
qp_matix = pickle.load(open('qp_matrix.pkl', 'rb')) 

c_phase_gate = Qobj(cp_matix,
              dims = [[2,2], [2,2]])

q_phase_gate = Qobj(qp_matix,
              dims = [[2,2], [2,2]])
q_phase_gate = q_phase_gate/q_phase_gate.data[0,0]

simulation_option = {
    "simulation_time": t_width, # ns
    "simulation_step": 10000,
    "initial_state": [state_0101, state_0110, state_1001, state_1010] # Do multiple simulation
}

def run_sim(args):

    leakage_list = []; angle_cphase=[]; Cphase_melements= []
    Fd_avg = []; Fd_process = []; Fd_trace = []

    repeat_time = args.rt
    Q_factor = args.Q
    B_noise_std = args.Bnoise
    J_noise_std = (1e-3+args.jnoise)/np.sqrt(2)/np.pi/Q_factor
    J23_noise_std = (args.j23noise+1e-3)/np.sqrt(2)/np.pi/Q_factor
    if args.eigenInit == 0:
        simulation_option["initial_state"] = [state_dudu, state_duud, state_uddu, state_udud]

    Uphase_ideal_matrix = np.array([[1,0,0,0], 
                [0, 1, 0, 0],
                [0,0, 1, 0],
                [0,0,0, np.exp(1j*(np.pi))]])
    Uphase_ideal = Qobj(Uphase_ideal_matrix, dims= [[2,2], [2,2]])

    for rtime in tqdm(range(repeat_time)):
        pulse_sequence = [
            {
                'pulse_index': 0, # [m,n] represents $\Omega_{mn}$
                'type': "INT",
                'pulse_shape': "cosine",
                't_delay': 0, # unit in ns
                't_width': t_width, # unit in ns
                't_plateau':0 , # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': bisa_amp, # XY: Rabi freq; Z: biased frequency
                'q_index': [1,2],
            },
            {
                'pulse_index': 1,
                'type': "Z",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, B_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': 0,
            },
            {
                'pulse_index': 2,
                'type': "Z",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, B_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': 1,
            },
            {
                'pulse_index': 3,
                'type': "Z",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, B_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': 2,
            },
            {
                'pulse_index': 4,
                'type': "Z",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, B_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': 3,
            },
            {
                'pulse_index': 5,
                'type': "INT",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, J_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': [0,1],
            },
            {
                'pulse_index': 6,
                'type': "INT",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, J_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': [2,3],
            },
            {
                'pulse_index': 7,
                'type': "INT",
                "pulse_shape": "square",
                't_delay': 0, # unit in ns
                't_width': 0, # unit in ns
                't_plateau':simulation_option['simulation_time'], # unit in ns
                'freq': 0, # unit in GHz; Z pulse does not use it
                'phase': 0, # unit in rad; Z pulse does not use it
                'amplitude': np.random.normal(0, J23_noise_std), # XY: Rabi freq; Z: biased frequency
                'q_index': [1,2],
            }
        ]
        result_list, angle_list = _system_arb1.system_dynamics_mesolve(simulation_option, pulse_sequence)
        propa_list = _system_arb1.system_dynamics_propagator(simulation_option, pulse_sequence, do_progress_bar=None)
        lk_dummy = []
        for i in range(len(simulation_option["initial_state"])):
            lk = 1 - np.abs(((state_udud + state_dudu + state_duud + state_uddu ).dag()* result_list[i].states[-1]).data[0,0])
            lk_dummy.append(lk)
        leakage_list.append(np.average(lk_dummy))

        U = propa_list[-1] # get the Unitary 
        # Perform partial trace, tracing out the coupler degree of freedom
        slist = [
            state_dudu, state_duud, state_uddu, state_udud
        ]
        sdlist = [
            state_dudu.dag(), state_duud.dag(), state_uddu.dag(), state_udud.dag()
        ]
        dims = [len(sdlist), len(slist)];   umatrix = []
        for i in range(dims[0]):
            umatrix_row = []
            for j in range(dims[1]):
                umatrix_row.append(sdlist[i] * U * slist[j])
            umatrix.append(umatrix_row)

        Usim = Qobj(np.array(umatrix).reshape(dims[0],dims[1]), dims = [[int(np.sqrt(dims[0])), int(np.sqrt(dims[0]))], [int(np.sqrt(dims[1])), int(np.sqrt(dims[1]))]])
        Usim = Usim/Usim.data[0,0]
        Usim_p = Usim * q_phase_gate
        Uphase = Usim_p * c_phase_gate

        angle_cphase.append(np.angle(Uphase.data[3,3]))
        Cphase_melements.append([Uphase.data[jj,jj] for jj in range(4)])
        process_fd = np.abs(process_fidelity(Uphase, Uphase_ideal, normalize=True))
        avg_fd = np.abs(average_gate_fidelity(Uphase, Uphase_ideal))
        trace_fd = np.abs((1/4*(Uphase_ideal*Uphase.dag()).tr()))
        if process_fd <= 1: Fd_process.append(process_fd)
        if avg_fd <= 1: Fd_avg.append(avg_fd)
        if trace_fd <= 1: Fd_trace.append(trace_fd)
    data = [Fd_process, Fd_avg, Fd_trace, pulse_sequence, [J_noise_std, J23_noise_std], B_noise_std]
    trail = args.trail
    pickle.dump(data, open(f'../../../Data/CZ_{repeat_time}trail_35ns_{trail}.pkl', 'wb'))


def parse_arguments():
    """
    Passing example:
    python .\CZ_ST0_Ramsy.py --rt=128 --trail=8 --tau_list 500 600 30 --ep_list 10 10 1
    """
    parser = argparse.ArgumentParser()
    # parser.add_argument('--note', type=str, default="fuck"
    #                     , help = '')
    parser.add_argument('--rt', type=int, default=1000,
                        help='')
    parser.add_argument('--trail', type=str, default="-1",
                        help='')
    parser.add_argument('--j23noise', type=float, default=5.152e-3,
                        help='')
    parser.add_argument('--jnoise', type=float, default=0,
                        help='')
    parser.add_argument('--Q', type=int, default=20,
                        help='')
    parser.add_argument('--Bnoise', type=float, default=0.22e-3,
                        help='')
    parser.add_argument('--eigenInit', type=int, default=1,
                        help='')
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    # print(str(args.type))
    run_sim(args)