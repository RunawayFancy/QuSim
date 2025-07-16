# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np


PI = np.pi


class ChrgNoiseExchangeQD:
    
    def __init__(self, Jres: float = 0.1e-3*2*PI, lever_arm: float = 1):
        self.Jres = Jres
        self.lever_arm = lever_arm


    def Jexchange(self, v_barrier: np.ndarray):
        return self.Jres*np.exp(2*self.lever_arm*v_barrier) # unit in GHz


    def J2Vbarrier(self, J_exchange: np.ndarray):
        if np.min(J_exchange/self.Jres + 1)<=0:
            print(np.min(J_exchange/self.Jres + 1))
        return 1/(2*self.lever_arm) * np.log(J_exchange/self.Jres + 1)
            

    def dJdV(self, v_barrier: np.ndarray):
        return self.Jres*2*self.lever_arm*np.exp(2*self.lever_arm*v_barrier) # unit in GHz


    def tdbase_charge_noise(self, tlist: np.ndarray, waveform: np.ndarray):
        """
        Suggested method: 'sum'
        """
        Vbarrier_arr = self.J2Vbarrier(waveform+self.Jres)
        sensitivity_arr = self.dJdV(Vbarrier_arr)
        return sensitivity_arr


    def tranfofn_charge_noise(self, v_barrier: np.ndarray):
        """
        Suggested method: 'multiply'
        """
        return np.exp(2*self.lever_arm*v_barrier)