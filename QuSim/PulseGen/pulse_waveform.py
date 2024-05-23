# -*- coding: utf-8 -*-
"""
@author: Pan Shi, Jiheng Duan
"""
import numpy as np
import qusim.Instruments.tools as tools
from qusim.PulseGen.pulse_shape import pulse_shape_dic
from qusim.PulseGen.noise_gen import noise_gen

# cos in one period
# cos in half period
# cosh


class pulse_lib:

    def __init__(self, pulse):
        self.pulse = pulse
        self.pulse_index = str(pulse['pulse_index'])
        self.t_width = pulse["t_width"]
        self.t_plateau = pulse["t_plateau"]
        self.t_delay = pulse["t_delay"]
        self.ampli = pulse["amplitude"]
        self.pulse_shape = pulse["pulse_shape"]
        self.pulse_type = pulse["type"]
        self.freq = pulse.get("freq", 0)
        self.phase = pulse.get("phase", 0)
        self.noise_chan = pulse.get("noise", 0)

        self.t_g = self.t_width + self.t_plateau
        self.t_rising = self.t_width/2
    
    def get_pulse(self, simulation_option):
        tlist = np.linspace(0, simulation_option["simulation_time"], simulation_option["simulation_step"])
        delta_t = tlist[1] - tlist[0]
        pulse_shape = pulse_shape_dic[self.pulse_shape]
        drive_pulse = pulse_shape(tlist, self.pulse)
        if drive_pulse is None:
            raise ValueError("Invalid pulse shape: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
    
        if self.pulse_type in ['XY', 'Z', 'INT']:
            if 'DRAG_scale' in self.pulse or 'DRAG_delta' in self.pulse:
                # DRAG correction
                # Multiple DRAG, define a new function
                drag_scale_list, drag_delta_list, num_drag = self.DRAG_pulse_pend_list()
                for i in range(num_drag):
                    drive_pulse += self.DRAG(drive_pulse, drag_scale_list[i], drag_delta_list[i], delta_t)

                # Pulse detuning correction on
                # After all DRAG
            # else: 
            #     if 'DRAG_scale' in self.pulse or 'DRAG_delta' in self.pulse:
            #         raise ValueError("Missing DRAG scale or DRAG delta: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)

            if 'pulse_detuning' in self.pulse:
                center = (self.t_width + self.t_plateau)/2 + self.t_delay
                drive_pulse *= np.exp(1j*2*np.pi*self.pulse['pulse_detuning'] * (tlist - center))

            # Multiply the carrier part
            carrier = self.carrier(tlist, self.freq, self.phase)
            drive_pulse *= carrier
            drive_pulse = np.real(drive_pulse)
        else: ValueError("Invalid pulse type: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
        # Add noise
        if self.noise_chan != 0:
            drive_pulse = self.add_noise(np.real(drive_pulse), simulation_option)
        if 'offset' in self.pulse:
            drive_pulse += self.pulse["offset"]
        return drive_pulse

    def carrier(self, tlist, freq, phase = 0):
        return np.exp(-1j * (2 * np.pi * freq * tlist + phase))
    
    def DRAG(self, drive_pulse, drag_scale, drag_delta, delta_t):
        # if np.abs(drag_delta) > 1: raise ValueError('DRAG delta value is too small = {}'.format(drag_delta))
        return -1j*tools.grad(drive_pulse)/delta_t * drag_scale/(drag_delta * 2 * np.pi) 

    def DRAG_pulse_pend_list(self):
        try:
            drag_scale = self.pulse['DRAG_scale']
        except KeyError:
            drag_scale = 1
        try:
            drag_delta = self.pulse['DRAG_delta']
        except KeyError:
            raise ValueError("Missing DRAG delta: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
        
        # drag_scale = self.pulse['DRAG_scale']
        # drag_delta = self.pulse['DRAG_delta']

        if isinstance(drag_scale, list):
            if isinstance(drag_delta, list):
                if len(drag_scale) != len(drag_delta): raise ValueError("Miss matching between DRAG scale and delta: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)

                num_drag = len(drag_scale)
                drag_scale_list = drag_scale
                drag_delta_list = drag_delta
            else: 
                num_drag = len(drag_scale)
                drag_scale_list = drag_scale
                drag_delta_list = [drag_delta for j in range(num_drag)]
        else:
            if drag_scale < 0 or drag_scale > 1: raise ValueError("Invalid DRAG scale. pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
            if isinstance(drag_delta, list):
                num_drag = len(drag_delta)
                drag_scale_list = [drag_scale for j in range(num_drag)]
                drag_delta_list = drag_delta
            else:
                num_drag = 1
                drag_scale_list = [drag_scale]
                drag_delta_list = [drag_delta]
        return drag_scale_list, drag_delta_list, num_drag
    
    def add_noise(self, wf, simulation_option):
        for config in self.noise_chan:
            noise = noise_gen(0, simulation_option["simulation_time"], simulation_option["simulation_step"], config)
            wf += np.real(noise)
        return wf
    # def get_pulse_old(self):
    #     if self.pulse_type == "XY":
    #         pulse_shape_mapping = {
    #             "cos": self.cos_waveform
    #         }
    #     elif self.pulse_type == 'Z':
    #         pulse_shape_mapping = {
    #             "sin": self.sin_waveform,
    #             "cosh": self.cosh_waveform
    #         }
    #     else: ValueError("Invalid pulse type: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
    #     drive_pulse = pulse_shape_mapping.get(self.pulse_shape)
    #     if drive_pulse is None:
    #         raise ValueError("Invalid pulse shape: pulse_index = " + self.pulse_index + ', pulse_shape = ' + self.pulse_shape)
    #     return drive_pulse

    # def DRAG(self):
    #     try:
    #         self.pulse['DRAG_scale']
    #     except KeyError:
    #         do_drag = 0
    #     if do_drag:
    #         def DRAG_waveform(t_relative, t_laggy): return - self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase + np.pi/2) * np.pi * self.DRAG_param * np.sin(np.pi * (t_relative + t_laggy)/self.t_rising)
    #     else:
    #         def DRAG_waveform(t_relative, t_laggy): return 0
    #     return DRAG_waveform
    
    # def cos_waveform(self, t, args):
    #     t_relative = t - self.t_delay
    #     if t_relative < 0: 
    #         return 0
    #     if 0 <= t_relative <= self.t_rising:
    #         return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase)* (1 - np.cos(np.pi * t_relative / self.t_rising)) + self.DRAG_waveform(t_relative, 0)
    #     if self.t_rising < t_relative < self.t_g - self.t_rising:
    #         return self.ampli * np.cos(self.freq * np.pi *  2 * t_relative)
    #     if self.t_g - self.t_rising <= t_relative  <= self.t_g:
    #         return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative + self.phase)* (1 - np.cos(np.pi * ( t_relative - self.t_g + 2 * self.t_rising)/self.t_rising)) + self.DRAG_waveform(t_relative, - self.t_g + 2 * self.t_rising)
    #     if self.t_g < t_relative: 
    #         return 0
        
    # def sin_waveform(self, t, args):
    #     t_relative = t - self.t_delay
    #     if t_relative < 0:
    #         return 0
    #     if 0 <= t_relative <= self.t_rising:
    #         return self.ampli * np.sin(np.pi/2 * t_relative/self.t_rising)
    #     if self.t_rising < t_relative <= self.t_g - self.t_rising:
    #         return self.ampli
    #     if self.t_g - self.t_rising <= t_relative <= self.t_g:
    #         return self.ampli * np.sin(np.pi/2 * (t_relative - self.t_g + 2 * self.t_rising)/self.t_rising)
    #     if self.t_g < t_relative:
    #         return 0
    
    # def cosh_waveform(self, t, args):
    #     t_relative = t - self.t_delay
    #     zero1, zero2 = np.log(2-np.sqrt(3)), np.log(2+np.sqrt(3))
    #     if t_relative < 0:
    #         return 0
    #     if 0 <= t_relative <= self.t_rising:
    #         x = zero1 * (t_relative/self.t_rising - 1)
    #         return self.ampli * ( - np.cosh(x) + 2)
    #     if self.t_rising < t_relative < self.t_g - self.t_rising:
    #         return self.ampli
    #     if self.t_g - self.t_rising <= t_relative <= self.t_g:
    #         x = zero2 * (t_relative - self.t_g + self.t_rising)/self.t_rising
    #         return self.ampli * ( - np.cosh(x) + 2)
    #     if self.t_g < t_relative:
    #         return 0