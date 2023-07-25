import numpy as np

# cos in one period
# cos in half period
# cosh

class pulse_lib:

    def __init__(self, pulse):
        self.t_width = pulse["t_width"]
        self.t_plateau = pulse["t_plateau"]
        self.t_g = self.t_width + self.t_plateau
        self.t_rising = self.t_width/2
        self.t_delay = pulse["t_delay"]
        self.ampli = pulse["amplitude"]
        self.freq = pulse["freq"]
        self.pulse_shape = pulse["pulse_shape"]
        self.pulse_type = pulse["type"]
    
    def get_pulse(self):
        if self.pulse_type == "XY":
            pulse_shape_mapping = {
                "cos": self.cos_waveform
            }
        else:
            pulse_shape_mapping = {
                "sin": self.sin_waveform
            }
        drive_pulse = pulse_shape_mapping.get(self.pulse_shape)
        if drive_pulse is None:
            raise ValueError("Invalid pulse shape: " + self.pulse_shape)
        return drive_pulse
    
    def cos_waveform(self, t, args):
        t_relative = t - self.t_delay
        if t_relative <= 0: 
            return 0
        if 0 < t_relative < self.t_rising:
            return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t)* (1 - np.cos(np.pi * t_relative / self.t_rising))
        if self.t_rising <= t_relative <= self.t_g - self.t_rising:
            return self.ampli * np.cos(self.freq * np.pi * 2 * t_relative)
        if self.t_g - self.t_rising < t_relative  < self.t_g:
            return self.ampli/2 * np.cos(self.freq * np.pi * 2 * t_relative)* (1 - np.cos(np.pi * ( t_relative - self.t_g + 2 * self.t_rising)/self.t_rising))
        if self.t_g <= t_relative: 
            return 0
        
    def sin_waveform(self, t, args):
        t_relative = t - self.t_delay
        if t_relative <= 0:
            return 0
        if 0 < t_relative < self.t_rising:
            return self.ampli * np.sin(np.pi/2 * t_relative/self.t_rising)
        if self.t_rising <= t_relative <= self.t_g - self.t_rising:
            return self.ampli
        if self.t_g - self.t_rising < t_relative <self.t_g:
            return self.ampli * np.sin(np.pi/2 * (t_relative - self.t_g + 2 * self.t_rising)/self.t_rising)
        if self.t_g < t_relative:
            return 0