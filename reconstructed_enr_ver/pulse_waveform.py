# -*- coding: utf-8 -*-
"""
@author: Pan Shi, Jiheng Duan
"""
import numpy as np

# cos in one period
# cos in half period
# cosh

class pulse_lib:

    def __init__(self, pulse):
        self.pulse_index = pulse['pulse_index']
        self.t_width = pulse["t_width"]
        self.t_plateau = pulse["t_plateau"]
        self.t_g = self.t_width + self.t_plateau
        self.t_rising = self.t_width/2
        self.t_delay = pulse["t_delay"]
        self.ampli = pulse["amplitude"]
        if "freq" in pulse.keys(): self.freq = pulse["freq"]
        if "phase" in pulse.keys(): self.phase = pulse["phase"]
        self.pulse_shape = pulse["pulse_shape"]
        self.pulse_type = pulse["type"]
        if 'DRAG_param' in pulse.keys(): self.DRAG_param = pulse['DRAG_param']
        else: self.DRAG_param = False
        self.DRAG_waveform = self.DRAG()
    
    def get_pulse(self):
        if self.pulse_type == "XY":
            pulse_shape_mapping = {
                "cos": self.cos_waveform
            }
        elif self.pulse_type == 'Z':
            pulse_shape_mapping = {
                "sin": self.sin_waveform,
                "cosh": self.cosh_waveform
            }
        else: ValueError("Invalid pulse type: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
        drive_pulse = pulse_shape_mapping.get(self.pulse_shape)
        if drive_pulse is None:
            raise ValueError("Invalid pulse shape: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
        return drive_pulse
    
    def DRAG(self):
        if self.DRAG_param:
            def DRAG_waveform(t_relative, t_laggy): return - self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase + np.pi/2) * np.pi * self.DRAG_param * np.sin(np.pi * (t_relative + t_laggy)/self.t_rising)
        else:
            def DRAG_waveform(t_relative, t_laggy): return 0
        return DRAG_waveform
    
    def cos_waveform(self, t, args):
        t_relative = t - self.t_delay
        if t_relative < 0: 
            return 0
        if 0 <= t_relative <= self.t_rising:
            return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase)* (1 - np.cos(np.pi * t_relative / self.t_rising)) + self.DRAG_waveform(t_relative, 0)
        if self.t_rising < t_relative < self.t_g - self.t_rising:
            return self.ampli * np.cos(self.freq * np.pi *  2 * t_relative)
        if self.t_g - self.t_rising <= t_relative  <= self.t_g:
            return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase)* (1 - np.cos(np.pi * ( t_relative - self.t_g + 2 * self.t_rising)/self.t_rising)) + self.DRAG_waveform(t_relative, - self.t_g + 2 * self.t_rising)
        if self.t_g < t_relative: 
            return 0
        
    def sin_waveform(self, t, args):
        t_relative = t - self.t_delay
        if t_relative < 0:
            return 0
        if 0 <= t_relative <= self.t_rising:
            return self.ampli * np.sin(np.pi/2 * t_relative/self.t_rising)
        if self.t_rising < t_relative <= self.t_g - self.t_rising:
            return self.ampli
        if self.t_g - self.t_rising <= t_relative <= self.t_g:
            return self.ampli * np.sin(np.pi/2 * (t_relative - self.t_g + 2 * self.t_rising)/self.t_rising)
        if self.t_g < t_relative:
            return 0
    
    def cosh_waveform(self, t, args):
        t_relative = t - self.t_delay
        zero1, zero2 = np.log(2-np.sqrt(3)), np.log(2+np.sqrt(3))
        if t_relative < 0:
            return 0
        if 0 <= t_relative <= self.t_rising:
            x = zero1 * (t_relative/self.t_rising - 1)
            return self.ampli * ( - np.cosh(x) + 2)
        if self.t_rising < t_relative < self.t_g - self.t_rising:
            return self.ampli
        if self.t_g - self.t_rising <= t_relative <= self.t_g:
            x = zero2 * (t_relative - self.t_g + self.t_rising)/self.t_rising
            return self.ampli * ( - np.cosh(x) + 2)
        if self.t_g < t_relative:
            return 0