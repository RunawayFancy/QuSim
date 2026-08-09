#%%

import sys
import copy
# This on is a path correcting code, which is used to find the path of qusim.
sys.path.append('../..')

from qusim import *
from qusim.data_plot.plot_tomo import hinton as hintonpt
from qusim.utils.noise_trafofn_tdbasefn import *
from qusim.utils.noise_debug import *

import numpy as np; PI = np.pi
import matplotlib.pyplot as plt

from tqdm import *
from time import *
from sympy import*
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import pandas as pd
import math as m
from typing import List
from qutip import*

import scipy as sp
import numpy as np



B_lst = [10.75, 10.6] #GHz
Jres = 0.1e-3
num_q = len(B_lst)
freq_list = [[-B_lst[i]/2 , B_lst[i]/2] for i in range(num_q)]
inter_list = [
    [{
        "v01": 0,
        "v00": 1,
        "v11": -1
    } for i in range(num_q)],
    [{
        "v01": 1,
        "v00": 0,
        "v11": 0
    } for i in range(num_q)],
    [{
        "v01": -1j,
        "v00": 0,
        "v11": 0
    } for i in range(num_q)]]
r_dic = {
    "r12": Jres
}
r = r2matrix(r_dic, freq_list) 
extra_list=None
driving_list=[{'W01':1j} for i in range(num_q)]
bias_list=[{"Z11":1, 'Z00':-1} for i in range(num_q)]
gamma_list = None
qbsystem = aqs.ArbQubitSys(freq_list, inter_list, r, extra_list, gamma_list, driving_list, bias_list)

# |\uparrow> = |1>; |\downarrow> = |0>;
state_01, E01, idx01 = qbsystem.get_eigenstates_energy((0,1))
state_10, E10, idx10 = qbsystem.get_eigenstates_energy((1,0))

state_00, E00, idx00 = qbsystem.get_eigenstates_energy((0,0))
state_11, E11, idx11 = qbsystem.get_eigenstates_energy((1,1))

state_ud = Qobj(np.array([0,1,0,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)
state_du = Qobj(np.array([0,0,1,0]), dims = [[2, 2], [1, 1]], shape = (4, 1), type = ket)
tau_larmor = 1/np.sqrt(np.abs(B_lst[0]-B_lst[1])**2 + 16*Jres**2)
cosT = np.abs(B_lst[0]-B_lst[1])*tau_larmor
sinT = 4*Jres*tau_larmor

from typing import List
def Xcompiler(num_of_p: int, params: dict, t_wait:float = 0) -> List[PulseConfig]:
    pp = PulseConfig(
        pulse_index=1,
        pulse_type='INT',
        pulse_shape=PulseShapeFn.COSINE,
        t_delay=0,
        t_width=0,
        t_plateau=0,
        qindex=[0,1],
        amplitude=0,
        frequency=0
    )
    pseq = []
    vz_x = params.get("vz", 0)
    vz_phi = vz_x/2
    
    t_moment_tri = params.get("t_moment", False)
    tau_larmor = params.get("tau_larmor", False)
    if t_moment_tri:
        t_wait += m.ceil(params["tw"]/tau_larmor)*tau_larmor-params["tw"]
    for _i in range(num_of_p):  
        pp.t_delay = _i * (params["tw"]+t_wait)
        pp.t_width = params["tw"]
        pp.amplitude = params["amp"]
        pp.frequency = params["freq"]
        pp.pulse_index = 2*_i+1
        pp.phase = vz_phi
        pseq.append(copy.deepcopy(pp))

        vz_phi += vz_x

        pp.t_delay = _i * (params["tw"]+t_wait)
        pp.t_width = params["tw"]
        pp.amplitude = params["amp"]
        pp.frequency = 0
        pp.pulse_index = 2*_i+2
        pp.phase = 0
        pseq.append(copy.deepcopy(pp))
    return pseq, (num_of_p * (params["tw"]+t_wait) - t_wait)

wd_larmor = 1/tau_larmor
tw_larmor = 6*tau_larmor
wd = wd_larmor
tw_larmor, wd_larmor

#%%

simopt = SimulationOption(init_state=[], simu_point=1000000, simu_time=100);
A_mu = 1e-6 # eV/sqrt(Hz)
# lever_arm_eff = 0.011 # 1/mV
lever_arm_eff = 0.1 # meV/mV
amp  = A_mu**2/(lever_arm_eff**2)

alpha = 1         # For 1/f noise

N_samples = simopt.simulation_point
t_start = 0  
t_end = simopt.simulation_time # unit in ns

f_min = 1e-3# unit Hz

time_series = np.linspace(t_start, t_end, N_samples)  # Convert ns to seconds
np.random.seed(1235456)
chrg_exchange = ChrgNoiseExchangeQD(Jres*2*PI, 0.011)
chrg_exchange_model = OneOverFNoiseConfig(
                                NoiseTimeConfig(simopt, tranfofn=chrg_exchange.tranfofn_charge_noise),
                                methods='multiply', 
                                lfreq=f_min, 
                                alpha=alpha, 
                                scale=amp,
                                corr_id = f'int_charge_{1}'
)



chrg_exchange_model2 = OneOverFNoiseConfig(
                                NoiseTimeConfig(simopt, tranfofn=chrg_exchange.tranfofn_charge_noise),
                                methods='multiply', 
                                lfreq=f_min, 
                                alpha=alpha, 
                                scale=amp,
                                corr_id = f'int_charge_{1}'
)

#%%

channel_noise = []
channel_noise.append(
    [
        ("INT", [1, 2]), 
        [chrg_exchange_model]
    ]
)
channel_noise.append(
    [
        ("INT", [2,3]), 
        [chrg_exchange_model2]
    ]
)

#%%

series = build_channel_noise_series(channel_noise,simopt)

# %%
noise1 = series[str(("INT", [1, 2]))]["raw"]
noise2 = series[str(("INT", [2, 3]))]["raw"]
frequencies, psd = welch(noise1, fs=chrg_exchange_model.sampling_freq, noverlap=0, nperseg=10000)
frequencies2, psd2 = welch(noise2, fs=chrg_exchange_model.sampling_freq, noverlap=0, nperseg=10000)

#%%
figure, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(time_series, noise1, label='Noise 1')
ax[1].plot(time_series, noise2, label='Noise 2')
ax[0].set_xlabel('Time (ns)')
ax[1].set_xlabel('Time (ns)')

# %%
f,corr = sp.signal.csd(noise1, noise2,fs=chrg_exchange_model.sampling_freq,noverlap=0,nfft=10000)

plt.plot(f,np.abs(corr)/np.sqrt(psd*psd2))
# %%
