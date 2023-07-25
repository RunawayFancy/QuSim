# -*- coding: utf-8 -*-
"""
@author: Ji Chu
"""
import math
import qutip as qt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np

def Qflatten(Q):
    return qt.Qobj(Q.full())

def generateBasicOperator(nTrunc):
    # generate basic operators. matrix truncated at nTrunc 
    I = qt.qeye(nTrunc)
    a = qt.destroy(nTrunc)
    x = a + a.dag()
    p = -1j*(a - a.dag())
    aa = a.dag() * a
    aaaa = a.dag() * a.dag() * a * a
    return {'I':I, 'a':a, 'x':x, 'p':p, 'aa':aa, 'aaaa':aaaa}

def Trunc_to_lowest_two_level(Q_state,ntrunc):
    # n qubit m levels ,|a0a1a2....an> ,totol dimension m^n=a0*m^(n-1)+a1*m^(n-2)+....ak*m^(n-k-1),a_i<m
    #
    Q_trunc_num=[]
    Q_state_matrix=Q_state.full()
    
    for num in range(len(Q_state_matrix)):
        i=num
        while (i>=ntrunc-1):
            if i % ntrunc >= 2:
                Q_trunc_num.append(num)
                break
            i=i // ntrunc

    for position_num in Q_trunc_num[::-1]:
        if Q_state.type=='ket':
            Q_state_matrix=np.delete(Q_state_matrix,position_num,0)         
        if Q_state.type=='bra':
            Q_state_matrix=np.delete(Q_state_matrix,position_num,1)          
        if Q_state.type=='oper':
            Q_state_matrix=np.delete(Q_state_matrix,position_num,1)  
            Q_state_matrix=np.delete(Q_state_matrix,position_num,0)              
    return qt.Qobj(Q_state_matrix)


def Solve_factor(rho,basis):
    rho_16_16=[]
    for vector in basis:
        rho_16_16.append(vector.full().flatten())
    rho_16_16_T=np.mat(np.array(rho_16_16).T)
    rho_T=np.mat(rho.full().flatten()).T
    basis_factor=rho_16_16_T.T*rho_T
    return np.asarray(basis_factor)

def calcaulate_belta_256_256(init_state_list,pauli_basis_list):
    belta_mnjk=np.zeros([256,256])*1j
    for m in range(16):
        for n in range(16):
            for j in range(16):
                epsilon_rho=pauli_basis_list[m]*init_state_list[j]*pauli_basis_list[n].dag()
                basis_factor=Solve_factor(epsilon_rho,init_state_list)
                for k in range(16):
                    belta_mnjk[j*16+k,m*16+n]=basis_factor[k][0]
#                     belta_mnjk[m*16+n,j*16+k]=basis_factor[k][0]
    return belta_mnjk

def calcaulate_lambda_16_16(measured_state_list,init_state_basis):
    lambda_factor=np.array([])
    for measured_state in measured_state_list:
        lambda_factor=np.append(lambda_factor,Solve_factor(measured_state,init_state_basis).T[0])
    return lambda_factor


List_sPauli = ['I','X','Y','Z']
List_mPauli = [qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()]
process_pauli_dict_16 =[]
for k1 in range(4):
    for k2 in range(4):
        process_pauli_dict_16.append(Qflatten(qt.tensor(List_mPauli[k1], List_mPauli[k2])))


class Two_qubit_process_tomo():

    def __init__(self,init_states_list,meas_state_list):
        self.pauli_basis_list=process_pauli_dict_16
        self.init_state_list=init_states_list
        self.measured_state_list=meas_state_list
        self.calcaulate_kappa_16_16()
    
    def calcaulate_kappa_16_16(self):
        self.lambda_factor_m=np.mat(calcaulate_lambda_16_16(self.measured_state_list,self.init_state_list)).T
        self.belta_m=np.mat(calcaulate_belta_256_256(self.init_state_list,self.pauli_basis_list))
        self.kappa_m=self.belta_m.I*self.lambda_factor_m
        self.kappa_list=np.asarray(self.kappa_m).T[0]

    def plot_bar3D(self):
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')
        _x = np.arange(16)
        _y = np.arange(16)
        _xx, _yy = np.meshgrid(_x, _y)
        x, y = _xx.ravel(), _yy.ravel()

        real_part =np.real(self.kappa_list)
        imag_part =np.imag(self.kappa_list)

        bottom1,top1,bottom2,top2 = np.zeros_like(real_part),np.zeros_like(real_part),np.zeros_like(real_part),np.zeros_like(real_part)

        for i in range(len(real_part)):
            top1[i]=max(real_part[i],0)
            bottom1[i]=min(real_part[i],0)
            top2[i]=max(imag_part[i],0)
            bottom2[i]=min(imag_part[i],0)

        width = depth = 0.85
        ax1.bar3d(x, y, bottom1, width, depth, top1, shade=True)
        ax2.bar3d(x, y, bottom2, width, depth, top2, shade=True)
        ax1.set_zlim([-1,1])
        ax2.set_zlim([-1,1])
        plt.show()






def calcaulate_belta_16_16(init_state_list,pauli_basis_list):
    belta_mnjk=np.zeros([16,16])*1j
    for m in range(4):
        for n in range(4):
            for j in range(4):
                epsilon_rho=pauli_basis_list[m]*init_state_list[j]*pauli_basis_list[n].dag()
                basis_factor=Solve_factor(epsilon_rho,init_state_list)
                for k in range(4):
                    belta_mnjk[j*4+k,m*4+n]=basis_factor[k][0]
#                     belta_mnjk[m*16+n,j*16+k]=basis_factor[k][0]
    return belta_mnjk

def calcaulate_lambda_4_4(measured_state_list,init_state_basis):
    lambda_factor=np.array([])
    for measured_state in measured_state_list:
        lambda_factor=np.append(lambda_factor,Solve_factor(measured_state,init_state_basis).T[0])
    return lambda_factor


List_sPauli = ['I','X','Y','Z']
List_mPauli = [qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()]
process_pauli_dict_4 = List_mPauli

class Single_qubit_process_tomo():

    def __init__(self,init_states_list,meas_state_list):
        self.pauli_basis_list=process_pauli_dict_4
        self.init_state_list=init_states_list
        self.measured_state_list=meas_state_list
        self.calcaulate_kappa_4_4()
    
    def calcaulate_kappa_4_4(self):
        self.lambda_factor_m=np.mat(calcaulate_lambda_4_4(self.measured_state_list,self.init_state_list)).T
        self.belta_m=np.mat(calcaulate_belta_16_16(self.init_state_list,self.pauli_basis_list))
        self.kappa_m=self.belta_m.I*self.lambda_factor_m
        self.kappa_list=np.asarray(self.kappa_m).T[0]

    def plot_bar3D(self):
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')
        _x = np.arange(16)
        _y = np.arange(16)
        _xx, _yy = np.meshgrid(_x, _y)
        x, y = _xx.ravel(), _yy.ravel()

        real_part =np.real(self.kappa_list)
        imag_part =np.imag(self.kappa_list)

        bottom1,top1,bottom2,top2 = np.zeros_like(real_part),np.zeros_like(real_part),np.zeros_like(real_part),np.zeros_like(real_part)

        for i in range(len(real_part)):
            top1[i]=max(real_part[i],0)
            bottom1[i]=min(real_part[i],0)
            top2[i]=max(imag_part[i],0)
            bottom2[i]=min(imag_part[i],0)

        width = depth = 0.85
        ax1.bar3d(x, y, bottom1, width, depth, top1, shade=True)
        ax2.bar3d(x, y, bottom2, width, depth, top2, shade=True)
        ax1.set_zlim([-1,1])
        ax2.set_zlim([-1,1])
        plt.show()

