from qusim import *
from qusim.Utils.noise_trafofn_tdbasefn import *
path = 'Data/ST0_coupler/'
qsv = QSave(path)

import numpy as np; PI = np.pi
import matplotlib.pyplot as plt

from tqdm import *
from time import *
from sympy import*
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import copy
import pickle
import math as m

from qutip import*

import argparse
from qusim.PulseGen.noise_gen import clear_corr_noise_cache

## Define system

B_lst = [10.75, 10.6, 10.54, 10.48, 10.35] #GHz
Jres_intra = 0.5e-3
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
extra_list=None
driving_list=[{'W01':1j} for i in range(num_q)]
bias_list=[{"Z11":1, 'Z00':-1} for i in range(num_q)]
gamma_list = None


def Xcompiler(num_of_p: int, params: dict, t_wait:float = 0) -> List[PulseConfig]:

    pseq = []
    vz_x = params.get("vz", 0)
    vz_phi = vz_x/2
    qindex  = params.get("qindex", [0,1])
    
    pp = PulseConfig(
        pulse_index=1,
        pulse_type='INT',
        pulse_shape=PulseShapeFn.COSINE,
        t_delay=0,
        t_width=0,
        t_plateau=0,
        qindex=qindex,
        amplitude=0,
        frequency=0
    )
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


def averaging_noisy_unitary(propagator_sublist: List[Qobj], initial_state: Qobj):
    N = len(propagator_sublist)
    empt_lst = []
    for _i, U in enumerate(propagator_sublist):
        empt_lst.append(U * initial_state * U.dag())
    return sum(empt_lst)/N

def transfer_basis(Mint: Qobj, tbasis:List[Qobj] ):
    Msub = np.zeros([len(tbasis), len(tbasis)],dtype='complex128')
    for _i, srow in enumerate(tbasis):
        for _j, scol in enumerate(tbasis):
            Msub[_i,_j] = Mint.matrix_element(srow.dag(), scol)
            # Msub[_i,_j] = (srow.dag() * Mint * scol).data[0,0]

    M_qobj = Qobj( Msub, dims = [[len(tbasis)], [len(tbasis)]], shape=[len(tbasis), len(tbasis)])
    return M_qobj

def transfer_basis_propagator_sublist(propagator_sublist, tbasis:List[Qobj], drive_correction, transfer_dim = None):
    res_lst = []
    for U in propagator_sublist:
        transfered_U = transfer_basis(U, tbasis)
        if transfer_dim is not None:
            transfered_U.dims = transfer_dim
        res_lst.append( transfered_U)
    return drive_correction *res_lst

def single_q_phase_correction(tw, qbsystem):
    corre_phase_matrix = (1j*tw*qbsystem.H).expm()

    return corre_phase_matrix

load_group_total = [
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.1MHz_Jres_inter_685.pkl", # loadID = 0
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.2MHz_Jres_inter_662.pkl", # loadID = 1
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.3MHz_Jres_inter_664.pkl", # loadID = 2
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.4MHz_Jres_inter_666.pkl", # loadID = 3
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.5MHz_Jres_inter_668.pkl", # loadID = 4
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.6MHz_Jres_inter_669.pkl", # loadID = 5
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.7MHz_Jres_inter_670.pkl", # loadID = 6
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.8MHz_Jres_inter_671.pkl", # loadID = 7
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_0.9MHz_Jres_inter_672.pkl", # loadID = 8
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_1.0MHz_Jres_inter_673.pkl", # loadID = 9
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_1.5MHz_Jres_inter_674.pkl", # loadID = 10
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_2.0MHz_Jres_inter_675.pkl", # loadID = 11
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_2.5MHz_Jres_inter_676.pkl", # loadID = 12
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_3.0MHz_Jres_inter_677.pkl", # loadID = 13
    "SpinCupler_Q1_2D_tw_ngate_repeat_gate_wo_tmoment_ctw_4.0MHz_Jres_inter_678.pkl", # loadID = 14
    ]

Jres_inter_lst_total = np.array([
                            0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3,
                            0.5e-3, 0.6e-3, 0.7e-3, 0.8e-3, 0.9e-3, 
                            1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 4e-3
                          ])
Jres_inter_scan_bound = [None, None, len(Jres_inter_lst_total)]

print(len(load_group_total), len(Jres_inter_lst_total))
assert len(load_group_total) == len(Jres_inter_lst_total)

# chrg_noise_trigger = False
# nrep = 1

do_save = True

tw_cd = 10 # unit in tau_larmor
qidx = [0,1]
ngate = 1
hyperfine_noise_trigger = True

sigma_bz_scan_bound = [0 , 0.2e-3, 10]
sigma_bz_scan_rng = np.linspace(*sigma_bz_scan_bound)

drive_oper = {"[0,1]": tensor(sigmaz(), qeye(2)), "[2,3]": tensor(qeye(2), sigmaz())}

# noise config
A_mu = 1e-6 # eV/sqrt(Hz)
lever_arm_eff = 0.1 # 1 mev/mV
amp  = A_mu**2/(lever_arm_eff**2)
alpha = 0.9         # For 1/f noise
chrg_exchange = ChrgNoiseExchangeQD(Jres_intra*2*PI, 0.011)

# sigma_bz = 0.06e-3#0.22e-3 # unit in GHz
bz_rescale_amp = 2*PI # unit in GHz
f_min = 1e-3# unit Hz

def run_simulation(args):
    
    load_group = [load_group_total[args.id]]
    Jres_inter_verify = Jres_inter_lst_total[args.id]
    chrg_noise_trigger = args.chrgtriger
    nrep = args.nrep
    do_scan = args.doscan
    # cntrl_handle = args.cntrl

    print(f"Jres_inter: {np.round(Jres_inter_verify*1e3,4)}MHz")
    print(f"Chrg noise trigger: {chrg_noise_trigger}")
    print(f"Hyperfine noise trigger: {hyperfine_noise_trigger}")
    print(f"nrep: {nrep}")
    print(f"load_group: {load_group[0]}")
    # print(f"control handle: {cntrl_handle}")

    # cntrl_handle_rng = [cntrl_handle]

    if do_scan:
        d_lst,fn_lst=qsv.load(load_group)
        propgator_lst = []
        for _i in tqdm(range(len(load_group))):
            d=d_lst[_i]
            fn=fn_lst[_i]
            ngate_scan_rng, tw_scan_rng, pop_lst = d['data']
            ngate_scan_bound, tw_scan_bound = d['rng']
            optimal_amp_list,vz_arr = d['other']

            Jres_inter = Jres_inter_verify
            tau_larmor = tw_scan_rng[0]/tw_scan_bound[0]
            try:
                opt_idx = np.where(np.linspace(*tw_scan_bound)==tw_cd)[0][0]
            except:
                raise ValueError("tw_cd not in the scan range")
            
            wd = 1/tau_larmor
            j0 = optimal_amp_list[opt_idx]
            tw = tw_scan_rng[opt_idx]
            vz = -vz_arr[opt_idx]
            param = {"amp": j0, "freq": wd, "tw": tw, 'vz': vz, "qindex": qidx}
            pseq, t_total = Xcompiler(ngate,param)
            npts = np.max([int(30*t_total) , 1000])
            simopt = SimulationOption(simu_time=t_total, simu_point=npts, init_state=[])

            chan_noise = []

            if chrg_noise_trigger:
                chrg_exchange_model = OneOverFNoiseConfig(
                    NoiseTimeConfig(simopt, tranfofn=chrg_exchange.tranfofn_charge_noise), 
                    methods='multiply', 
                    lfreq=f_min, 
                    alpha=alpha, 
                    scale=amp,
                    corr_id = 'int_charge'
                )
                chan_noise.append([("INT", [0,1]), [chrg_exchange_model]])
                chan_noise.append([("INT", [3,4]), [chrg_exchange_model]])

                chrg_exchange_inter = ChrgNoiseExchangeQD(Jres_inter*2*PI, 0.011)
                chrg_exchange_inter_model = OneOverFNoiseConfig(
                    NoiseTimeConfig(simopt, tranfofn=chrg_exchange_inter.tranfofn_charge_noise), 
                    methods='multiply', 
                    lfreq=f_min, 
                    alpha=alpha, 
                    scale=amp,
                    corr_id = 'int_charge'
                )
                chan_noise.append([("INT", [1,2]), [chrg_exchange_inter_model]])
                chan_noise.append([("INT", [2,3]), [chrg_exchange_inter_model]])
            else:
                chrg_exchange_model = None

            r_dic = {
                "r12": Jres_intra,
                "r45": Jres_intra,
                "r23": Jres_inter,
                "r34": Jres_inter
            }
            r = r2matrix(r_dic, freq_list) 
            qbsystem = aqs.ArbQubitSys(freq_list, inter_list, r, extra_list, gamma_list, driving_list, bias_list)

            s_00, E_00, idx_00 = qbsystem.get_eigenstates_energy((0,1,0,0,1)) # 00
            s_01, E_01, idx_01 = qbsystem.get_eigenstates_energy((0,1,0,1,0)) # 01
            s_10, E_10, idx_10 = qbsystem.get_eigenstates_energy((1,0,0,0,1)) # 10
            s_11, E_11, idx_11 = qbsystem.get_eigenstates_energy((1,0,0,1,0)) # 11

            if np.abs(min(s_10)[0][0]) > 0.9:
                s_10 = -s_10
            if np.abs(min(s_01)[0][0]) > 0.9:
                s_01 = -s_01
            if np.abs(min(s_11)[0][0]) > 0.9:
                s_11 = -s_11
            if np.abs(min(s_00)[0][0]) > 0.9:
                s_00 = -s_00
            assert np.abs(E_10-E_00 - wd)<1e-5


            propgator_lst.append([])

            for _j, sigma_bz in tqdm(enumerate(sigma_bz_scan_rng), leave=False):
                assert hyperfine_noise_trigger
                chan_noise_scan = copy.deepcopy(chan_noise)
                hyperfin_model = GaussianNoiseConfig(
                    NoiseTimeConfig(simopt), 
                    mean=0, 
                    std=sigma_bz, 
                    amp=bz_rescale_amp,
                    corr_id = 'hyperfine_z'
                    )# 2pi/2 means 1/2 \omega 2 pi sigma_z
                chan_noise_scan.append([("Z", 0), [hyperfin_model]])
                chan_noise_scan.append([("Z", 1), [hyperfin_model]])
                chan_noise_scan.append([("Z", 2), [hyperfin_model]])
                chan_noise_scan.append([("Z", 3), [hyperfin_model]])
                chan_noise_scan.append([("Z", 4), [hyperfin_model]])

                dummy_prop_lst = []

                # Pre-generate per-rep seeds so each noise realization is distinct
                rep_seeds = np.random.default_rng().integers(0, 2**32 - 1, size=nrep, dtype=np.uint32)

                for _k in range(nrep):
                    # Reset correlated cache and RNG so each rep uses a fresh correlated draw
                    clear_corr_noise_cache()
                    np.random.seed(int(rep_seeds[_k]))
                    U = qbsystem.system_dynamics_propagator(pseq, simopt, channel_noise=chan_noise_scan)[-1]
                    dummy_prop_lst.append(U * single_q_phase_correction(tw, qbsystem))
                qdidx = param["qindex"]
                drive_correction = (-1j*vz*drive_oper[f"[{qdidx[0]},{qdidx[1]}]"]/2).expm()
                U_reduced_lst = transfer_basis_propagator_sublist(dummy_prop_lst, [s_00, s_01, s_10, s_11], drive_correction, [[2,2],[2,2]])
                propgator_lst[_i].append(U_reduced_lst)


        if chrg_noise_trigger:
            chrg_word = 'w'
        else:
            chrg_word = 'wo'
        if hyperfine_noise_trigger:
            hyperfine_word = 'w'
        else:
            hyperfine_word = 'wo'

        if do_save:
            data = {"data": [sigma_bz_scan_rng, Jres_inter_lst_total, propgator_lst], "rng": [sigma_bz_scan_bound, Jres_inter_scan_bound, ngate, nrep, tw_cd, param]}

            fn = qsv.save(f"SpinCoupler_2D_SQ_crosstalk_{chrg_word}_CorrChrgNoise_{hyperfine_word}_CorrHyperNoise", data =data)



def parse_arguments():
    """
    Passing example:
    python .\CZ_ST0_Ramsy.py --doscan=0 --nrep=512 --id=0 --chrgtriger=0
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--nrep', type=int, default=512,
                        help='')
    parser.add_argument('--id', type=int, default=0,
                        help='')
    parser.add_argument('--chrgtriger', type=int, default=1,
                        help='')
    parser.add_argument('--doscan', type=int, default=0,
                        help='')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    # print(args.rt)
    # print(type(args.ep_list))
    # sys.exit(0)
    run_simulation(args)
