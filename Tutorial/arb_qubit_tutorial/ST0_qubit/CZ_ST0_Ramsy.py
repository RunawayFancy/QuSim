import sys
# This on is a path correcting code, which is used to find the path of qusim.
sys.path.append('../../..')
import qusim.DataPlot.plot_lib as pl
import qusim.Instruments.tools as tools
import qusim.Instruments.stochastic_gen as sg
import qusim.System.arb_qubit_system as aqs
import qusim.Instruments.angle as  ang
from qusim.Instruments.angle import get_angle

import numpy as np
import json
import matplotlib.pyplot as plt

from tqdm import *
from time import *
from sympy import*
import copy
import pickle
from qutip import*

import argparse

def var_name2str(variable):
    def nested_var_name2str(value):
        if isinstance(value, list):
            return [nested_var_name2str(item) for item in value]
        else:
            for name in globals():
                if eval(name) == value:
                    return name
            return str(value)
    return nested_var_name2str(variable)


freq_list = [ # GHz
    [0, 10],
    [0, 10.25]
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
    }]
]

r_dic = {
    "r12": 0e-3
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
    }
]

## with pulse type Z
bias_list = [ # Default, Z00 = 0
    {
        "Z11": 1
    },
    {
        "Z11": 1
    }
]

gamma_list = [{"z":1/50}, {"z": 1/50}]
gamma_list = None

_system_arb1 = aqs.arb_qubit_system(freq_list, inter_list, extra_list=extra_list, r =r, driving_list=driving_list, bias_list=bias_list, gamma_list=gamma_list)

# System initial state & eigenstate
state_00, E_00, _= _system_arb1.get_eigenstates_energy((0,0))
state_10, E_10, _= _system_arb1.get_eigenstates_energy((1,0))
state_01, E_01, _= _system_arb1.get_eigenstates_energy((0,1))
state_11, E_11, _= _system_arb1.get_eigenstates_energy((1,1))

# spin basis
state_uu = state_11
state_dd = state_00
state_ud = Qobj(np.array([0,0,1,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)
state_du = Qobj(np.array([0,1,0,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)

# ST0
state_S = (state_ud-state_du)/np.sqrt(2)
state_T0 = (state_ud+state_du)/np.sqrt(2)

init_state_dic ={
    "S": state_S,
    "T0": state_T0,
    "ud": state_ud,
    "du": state_du
}

def run_simulation(args):
    repeat_time = args.rt
    tau_list = np.linspace(*args.tau_list)
    # ep_list= np.linspace(100,0, 10)
    ep_list= np.linspace(*args.ep_list)
    j23_list = 10*np.exp(-ep_list/10)
    half_pi_gt = np.sqrt(2); bisa_amp = 61.5e-3
    noise_var = args.noise
    ramsy_pop = []
    noise_type = args.type
    init_state = init_state_dic[args.inistate]
    if noise_type == "Z":
        q_indeces = 1
    else:
        q_indeces = [0,1]

    for tau_d in tqdm(tau_list):
        total_t = tau_d + 2*half_pi_gt
        simulation_option = {
            "simulation_time": total_t, # ns
            "simulation_step": 10000,
            "initial_state": [init_state] # Do multiple simulation
        } 
        initial01_track01 = np.zeros(shape=(repeat_time,), dtype=np.float32)
        for tt in range(repeat_time):

            pulse_sequence = [
                {
                    'pulse_index': 0, # [m,n] represents $\Omega_{mn}$
                    'type': "INT",
                    'pulse_shape': "square",
                    't_delay': 0, # unit in ns
                    't_width': half_pi_gt, # unit in ns
                    't_plateau':0, # unit in ns
                    'freq': 0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': bisa_amp+1e-3,#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                    'q_index': [0,1]
                },
                {
                    'pulse_index': 1, # [m,n] represents $\Omega_{mn}$
                    'type': "INT",
                    'pulse_shape': "square",
                    't_delay': tau_d+half_pi_gt, # unit in ns
                    't_width': half_pi_gt, # unit in ns
                    't_plateau':0, # unit in ns
                    'freq': 0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': bisa_amp+1e-3,#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                    'q_index': [0,1]
                },
                { # Base line J23
                    'pulse_index': 4,
                    'type': 'Z',
                    "pulse_shape": "square",
                    't_delay': half_pi_gt, # unit in ns
                    't_width': 0, # unit in ns
                    't_plateau':tau_d, # unit in ns
                    'freq': 0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': np.random.normal(0, noise_var),#j23_d*1e-3, # XY: Rabi freq; Z: biased frequency
                    'q_index': 0
                },
                { # Base line J23
                    'pulse_index': 5,
                    'type': 'Z',
                    "pulse_shape": "square",
                    't_delay': half_pi_gt, # unit in ns
                    't_width': 0, # unit in ns
                    't_plateau':tau_d, # unit in ns
                    'freq': 0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': np.random.normal(0, noise_var),#j23_d*1e-3, # XY: Rabi freq; Z: biased frequency
                    'q_index': 1
                }
            ]

            result_list, _ = _system_arb1.system_dynamics_mesolve(simulation_option, pulse_sequence)
            initial01_track01[tt] = np.abs(((result_list[0].states[-1]).dag()*init_state)[0][0][0])**2
            # initial01_track01.append(expect(result_list[0].states[-1], ket2dm(state_ud)))
        avg_pop = np.average(initial01_track01)
        ramsy_pop.append(avg_pop)
    data = [ramsy_pop, tau_list, j23_list, ep_list, [half_pi_gt, repeat_time]]

    trail = args.trail
    pickle.dump(data, open(f'../../../Data/Qfactor_ramsy_{noise_type}noise_{trail}.pkl', 'wb'))


def parse_arguments():
    """
    Passing example:
    python .\CZ_ST0_Ramsy.py --rt=128 --trail=8 --tau_list 500 600 30 --ep_list 10 10 1
    """
    parser = argparse.ArgumentParser()
    # parser.add_argument('--note', type=str, default="fuck"
    #                     , help = '')
    parser.add_argument('--rt', type=int, default=1,
                        help='')
    parser.add_argument('--trail', type=int, default=-1,
                        help='')
    parser.add_argument('--noise', type=float, default=0.22e-3,
                        help='')
    parser.add_argument('--tau_list', nargs=3, type=int, default=[0, 100, 1000],
                        help='')
    parser.add_argument('--ep_list', nargs=3, type=int, default=[10, 10, 1],
                        help='')
    parser.add_argument('--type', type=str, default='Z',
                        help='')
    parser.add_argument('--inistate', type=str, default='ud',
                    help='ud, du, S, T0')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    # print(str(args.type))
    run_simulation(args)