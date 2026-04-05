#%%
import numpy as np
import matplotlib.pyplot as plt
import sys
from qutip import *
# # This on is a path correcting code, which is used to find the path of qusim.
sys.path.append('../..')

from qusim import *
from tqdm import *

qsv = QSave('../../Data/Test/')

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


#%% Define the system

N = None # Turn on maximum excitation level
# N = None # Turn off the maximum excitation level

w = [3] # Qubit frequency
q_dim = [2 for _ in range(len(w))] # Dimension of each qubit
alpha =[-0.3] # Anharmonicity
r = 0
gamma_list = [
    {
        "up": 0,
        "down": 0,
        "z": 0.05
    }
]

# Set up system class
_system = qs.TransmonSys(N, q_dim, w, alpha, r, gamma_list)
# Or
_system = qs.TransmonSys(N, q_dim, w, alpha)

# Get system Hamiltonian
_system.H
_system = qs.TransmonSys(N, q_dim, w, alpha, r, gamma_list)

state_0, E_0, index_0 = _system.get_eigenstates_energy((0,))
state_1, E_1, index_1 = _system.get_eigenstates_energy((1,))

#%% Define the pulse

t_half_pi = 20
A_half_pi = 1/4/t_half_pi
t_space = 10

wd_scan_rng = np.linspace(w[0]+0.025, w[0]+0.075, 64)
tplateau_scan_rng = np.linspace(0, 300, 64)

dt = 0.1 # ns
pop_1 = np.zeros([len(wd_scan_rng), len(tplateau_scan_rng)])

for _i, wd in tqdm(enumerate(wd_scan_rng)):
    for _j, tplateau in tqdm(enumerate(tplateau_scan_rng), disable=None):
        # print(f"wd: {wd:.3f}, tplateau: {tplateau:.1f}")
        t_total = 2*t_half_pi + t_space*2 + tplateau
        simopt = SimulationOption(
            simu_time=t_total,
            simu_point=0,   # ignored if dt is given
            dt=dt,
            init_state=[ket2dm(state_0)]
        )
        assert simopt.dt == dt, "The time step of the simulation option must be equal to the time step used for predistortion."
        pulse_sequence = [
            PulseConfig(
                pulse_index=0,
                pulse_type="XY",
                pulse_shape=PulseShapeFn.COSINE,
                t_delay=0,
                t_width=t_half_pi,
                t_plateau=0,
                frequency=wd,
                phase=0,
                amplitude=A_half_pi,
                qindex=0,
                # predistortion=[
                #     {"b": [2.48242606e-05, 4.96485211e-05, 2.48242606e-05],
                #     "a": [1.0, -1.9858581, 0.98595739]}
                # ]
            ),
            PulseConfig(
                pulse_index=1,
                pulse_type="Z",
                pulse_shape=PulseShapeFn.COSH,
                t_delay=t_space+t_half_pi,
                t_width=0,
                t_plateau=tplateau,
                frequency=0,
                phase=0,
                amplitude=0.05,
                qindex=0,
                epsilon=2,
                # predistortion=[
                #     {"b": [2.48242606e-05, 4.96485211e-05, 2.48242606e-05],
                #     "a": [1.0, -1.9858581, 0.98595739]}
                # ]
            ),
            PulseConfig(
                pulse_index=2,
                pulse_type="XY",
                pulse_shape=PulseShapeFn.COSINE,
                t_delay=t_space*2+t_half_pi+tplateau,
                t_width=t_half_pi,
                t_plateau=0,
                frequency=wd,
                phase=0,
                amplitude=A_half_pi,
                qindex=0,
                # predistortion=[
                #     {"b": [2.48242606e-05, 4.96485211e-05, 2.48242606e-05],
                #     "a": [1.0, -1.9858581, 0.98595739]}
                # ]
            ),
        ]
        # plot_pulse_sequence(pulse_sequence, simopt)
        result_list, _ = _system.system_dynamics_mesolve(pulse_sequence, simopt) 
        res = result_list[0]
        pop_1[_i,_j] = expect(res.states[-1], ket2dm(state_0))

#%%

data = {
    "wd_scan_rng": wd_scan_rng,
    "tplateau_scan_rng": tplateau_scan_rng,
    "pop_1": pop_1,
}

fn = qsv.save('2D_Ramsey_scan', data)

#%%
fig, ax = plt.subplots(1, figsize=(5,4))
pm = ax.pcolormesh(tplateau_scan_rng, wd_scan_rng, pop_1)
ax.set_ylabel(r"$\omega_d$ (GHz)", fontsize=12)
ax.set_xlabel("Plateau time (ns)", fontsize=12)
cb = plt.colorbar(pm, label=r"$\langle 0|\rho|0\rangle$", ax=ax)
plt.title(f"2D Ramsey {fn}")
plt.show()
# %%
