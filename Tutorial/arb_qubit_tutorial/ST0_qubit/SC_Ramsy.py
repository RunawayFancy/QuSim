import sys
# This on is a path correcting code, which is used to find the path of QuSim.
sys.path.append('../../..')
import QuSim.DataPlot.plot_lib as pl
import QuSim.Instruments.tools as tools
import QuSim.Instruments.stochastic_gen as sg
import QuSim.System.transmon_system as qs
import QuSim.Instruments.angle as  ang
from QuSim.Instruments.angle import get_angle

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


N = None # Maximum excitation level
w = [2] # Qubit frequency
q_dim = [2 for _ in range(len(w))]
alpha =[-0.3] # Anharmonicity

_system = qs.qubit_system(N, q_dim, w, alpha)

state_0, E_0, index_0 = _system.get_eigenstates_energy((0,))
state_1, E_1, index_1 = _system.get_eigenstates_energy((1,))

init_state_dic ={
    "0": state_0,
    "1": state_1,
}

def run_simulation(args):
    repeat_time = args.rt
    tau_list = np.linspace(*args.tau_list)
    half_pi_gt = 10; bisa_amp = 1/(2*half_pi_gt)
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
                    'type': "XY",
                    'pulse_shape': "cosine",
                    't_delay': 0, # unit in ns
                    't_width': half_pi_gt, # unit in ns
                    't_plateau':0, # unit in ns
                    'freq': E_1 - E_0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': bisa_amp,#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                    'q_index': 0
                },
                {
                    'pulse_index': 0, # [m,n] represents $\Omega_{mn}$
                    'type': "XY",
                    'pulse_shape': "cosine",
                    't_delay': half_pi_gt+tau_d, # unit in ns
                    't_width': half_pi_gt, # unit in ns
                    't_plateau':0, # unit in ns
                    'freq': E_1 - E_0, # unit in GHz; Z pulse does not use it
                    'phase': 0, # unit in rad; Z pulse does not use it
                    'amplitude': bisa_amp,#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                    'q_index': 0
                },
                # {
                #     'pulse_index': 0, # [m,n] represents $\Omega_{mn}$
                #     'type': "Z",
                #     'pulse_shape': "square",
                #     't_delay': 0, # unit in ns
                #     't_width': half_pi_gt, # unit in ns
                #     't_plateau':0, # unit in ns
                #     'freq': 0, # unit in GHz; Z pulse does not use it
                #     'phase': 0, # unit in rad; Z pulse does not use it
                #     'amplitude': 10.25 -freq_list[1][1],#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                #     'q_index': 1
                # },
                # {
                #     'pulse_index': 1, # [m,n] represents $\Omega_{mn}$
                #     'type': "Z",
                #     'pulse_shape': "square",
                #     't_delay': tau_d+half_pi_gt, # unit in ns
                #     't_width': half_pi_gt, # unit in ns
                #     't_plateau':0, # unit in ns
                #     'freq': 0, # unit in GHz; Z pulse does not use it
                #     'phase': 0, # unit in rad; Z pulse does not use it
                #     'amplitude': 10.25 -freq_list[1][1],#-3.5e-3, # bisa_amp, # XY: Rabi freq; Z: biased frequency
                #     'q_index': 1
                # },
                # { # Base line J23
                #     'pulse_index': 4,
                #     'type': 'Z',
                #     "pulse_shape": "square",
                #     't_delay': half_pi_gt, # unit in ns
                #     't_width': 0, # unit in ns
                #     't_plateau':tau_d, # unit in ns
                #     'freq': 0, # unit in GHz; Z pulse does not use it
                #     'phase': 0, # unit in rad; Z pulse does not use it
                #     'amplitude': np.random.normal(0, noise_var),#j23_d*1e-3, # XY: Rabi freq; Z: biased frequency
                #     'q_index': 0
                # },
                # { # Base line J23
                #     'pulse_index': 5,
                #     'type': 'Z',
                #     "pulse_shape": "square",
                #     't_delay': half_pi_gt, # unit in ns
                #     't_width': 0, # unit in ns
                #     't_plateau':tau_d, # unit in ns
                #     'freq': 0, # unit in GHz; Z pulse does not use it
                #     'phase': 0, # unit in rad; Z pulse does not use it
                #     'amplitude': np.random.normal(0, noise_var),#j23_d*1e-3, # XY: Rabi freq; Z: biased frequency
                #     'q_index': 1
                # }
            ]

            result_list, _ = _system.system_dynamics_mesolve(simulation_option, pulse_sequence)
            initial01_track01[tt] = np.abs( (init_state.dag() * result_list[0].states[-1]).data[0,0] )
            # initial01_track01.append(expect(result_list[0].states[-1], ket2dm(state_ud)))
        avg_pop = np.average(initial01_track01)
        ramsy_pop.append(avg_pop)
    data = [ramsy_pop, tau_list, 1, 1, [half_pi_gt, repeat_time]]

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