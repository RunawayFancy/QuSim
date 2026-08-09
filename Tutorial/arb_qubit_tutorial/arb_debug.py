#%% Import
import numpy as np
import matplotlib.pyplot as plt
import sys

import qusim.system.arb_qubit_system as aqs
# Some useful tools
import qusim.instruments.tools as tools
from qusim import *
from qutip import *

# Some intrinsic plotting function
import qusim.data_plot.plot_lib as pl

# Print the full output
np.set_printoptions(threshold=sys.maxsize)

# auto reload
%load_ext autoreload
%autoreload 2 

# variable name -> string
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

#%% Define Hamiltonian

freq_list = [
    [0, 6.3, 6],
    [0, 5.85, 5.58],
    [0, 6.2, 5.7]
]

inter_list = [
    {
        "v01": 1
    },
    {
        "v01": 1
    },
    {
        "v01": 1
    }
]

r_dic = {
    "r12": 0.05,
    "r13": -0.005,
    "r23": 0.05
}
r = tools.r2matrix(r_dic, freq_list) # Coupling strength

driving_list = [
    {
        "W01": 1,   # W01 represents \Omega_{01}, for the qubit with index 0
        "W12": np.sqrt(2)
    },
    {
        "W01": np.sqrt(2),
        "W12": np.sqrt(3)
    },
    {
        "W01": np.sqrt(3),
        "W12": np.sqrt(4)
    }
]

bias_list = [ # Default, Z00 = 0
    {
        "Z11": 1.1,
        "Z22": 2.1
    },
    {
        "Z11": 1.2,
        "Z22": 2.2
    },
    {
        "Z11": 1.3,
        "Z22": 2.3
    }
]

gamma_list = [
    {
        'up': 0.01,
        'down': 0.03,
        'z': 0.01
    },
    {
        'up': 0.01,
        'down': 0.1,
        'z': 0.01
    },
    {
        'up': 0.01,
        'down': 0.03,
        'z': 0.01
    }
]


# gamma_list = None

_system_arb = aqs.ArbQubitSys(freq_list, inter_list, r, gamma_list, driving_list, bias_list)

_system_arb.H

#%% Sys time evolution, use ket

freq_list = [
    [0, 6.3, 6.3*2 - 0.3],
    [0, 5.85, 5.85*2 - 0.3],
    [0, 6.2, 6.2*2 - 0.2]
]

inter_list = [
    {
        "v01": 1,
        "v12": np.sqrt(2)
    },
    {
        "v01": 1,
        "v12": np.sqrt(2)
    },
    {
        "v01": 1,
        "v12": np.sqrt(2)
    }
]

r_dic = {
    "r12": 0.05,
    "r13": -0.005,
    "r23": 0.05
}
r = tools.r2matrix(r_dic, freq_list) # Coupling strength

driving_list = [
    {
        "W01": 1,   # W01 represents \Omega_{01}, for the qubit with index 0
        "W12": np.sqrt(2)
    },
    {
        "W01": 1,
        "W12": np.sqrt(2)
    },
    {
        "W01": 1,
        "W12": np.sqrt(2)
    }
]

bias_list = [ # Default, Z00 = 0
    {
        "Z11": 1,
        "Z22": 2
    },
    {
        "Z11": 1,
        "Z22": 2
    },
    {
        "Z11": 1,
        "Z22": 2
    }
]

gamma_list = None

_system_arb = aqs.ArbQubitSys(freq_list, inter_list, r, gamma_list, driving_list, bias_list)

state_000, E_000, _ = _system_arb.get_eigenstates_energy((0,0,0))
state_001, E_001, _ = _system_arb.get_eigenstates_energy((0,0,1))
state_100, E_100, _ = _system_arb.get_eigenstates_energy((1,0,0))
state_101, E_101, _ = _system_arb.get_eigenstates_energy((1,0,1))
state_111, E_111, _ = _system_arb.get_eigenstates_energy((1,1,1))
state_200, E_200, _ = _system_arb.get_eigenstates_energy((2,0,0))

simopt = SimulationOption(
    simu_time=80,
    simu_point=10000,
    init_state=[state_000]
)

pseq = [
    PulseConfig(
        pulse_index=1,
        pulse_type = 'XY',
        pulse_shape=PulseShapeFn.COSINE,
        t_delay=10,
        t_width=50,
        t_plateau=0,
        frequency=(E_200 - E_000)/2,
        phase=0,
        amplitude=0.05,
        qindex=0,
    )
]

%matplotlib inline
## Notice that the all pulses'  amplitude are rescaled by a factor 1/1.2
pl.plot_pulse_sequence(pseq, simopt)

result_list, angle_list = _system_arb.system_dynamics_mesolve(pseq, simopt)


# state that you want to plot each simulation
interested_state = [
                    [state_200, state_100, state_000]
                    ]
interested_state_label = var_name2str(interested_state)
initial_state_label = var_name2str(simopt.initial_state)

%matplotlib inline
# plot state population evolution
pl.plot_population_evolution(_system_arb, result_list, simopt, interested_state, interested_state_label, initial_state_label)

#%% Use DM

state_000, E_000, index_000 = _system_arb.get_eigenstates_energy((0,0,0))
state_010, E_010, index_010 = _system_arb.get_eigenstates_energy((0,1,0))
state_001, E_001, index_001 = _system_arb.get_eigenstates_energy((0,0,1))
state_100, E_100, index_100 = _system_arb.get_eigenstates_energy((1,0,0))
state_101, E_101, index_101 = _system_arb.get_eigenstates_energy((1,0,1))
state_111, E_111, index_111 = _system_arb.get_eigenstates_energy((1,1,1))
state_200, E_200, index_200 = _system_arb.get_eigenstates_energy((2,0,0))
state_020, E_020, index_020 = _system_arb.get_eigenstates_energy((0,2,0))

simopt = SimulationOption(
    simu_time=80,
    simu_point=10000,
    init_state=[ket2dm(state_000)]
)

pseq = [
    PulseConfig(
        pulse_index=1,
        pulse_type = 'XY',
        pulse_shape=PulseShapeFn.COSINE,
        t_delay=30,
        t_width=30,
        t_plateau=0,
        frequency=(E_100 - E_000),
        phase=0,
        amplitude=0.2,
        qindex=0,
    )
]


%matplotlib inline
## Notice that the all pulses'  amplitude are rescaled by a factor 1/1.2
pl.plot_pulse_sequence(pseq, simopt)

result_list, angle_list = _system_arb.system_dynamics_mesolve(pseq, simopt)

t = simopt.tlist

# Select which result you want to see
result = result_list[0]

# Plot
%matplotlib inline

plt.plot(t,aqs.expect(result.states, state_000 * state_000.dag()), label=r'$\rho_{000,000}$');

plt.plot(t,aqs.expect(result.states, state_100 * state_100.dag()), label=r'$\rho_{100,100}$');

plt.ylabel(r"$\rho_{ij}$")
plt.xlabel("t (QuTiP 'seconds')")
plt.legend()
plt.show()
# %%
