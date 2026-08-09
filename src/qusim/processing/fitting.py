# -*- coding: utf-8 -*-
"""
@author: Zihao Wang
Original source: https://github.com/Zihao96/qtrlb/blob/main/qtrlb/processing/fitting.py
"""
# =============================================================================
# All the function in this script are supposed to be purely mathematical without
# considering the parameter or dictionary structure of measurement, so that it 
# could also be called for other purpose.
# 
# Example: Suppose we have 'data' as y axis and a known x axis called 't'.   
# 
#   fitmodel = ExpSinModel()
#   params = fitmodel.guess(data, t)
#   result = fitmodel.fit(data, params=params, x=t)
#
# The fitting result is stored in a dictionary result.best_values.
# =============================================================================

import numpy as np
from lmfit.models import Model, SineModel, QuadraticModel
from numpy import exp, sin
from scipy.signal import find_peaks
from qusim.utils.general_utils import make_it_list
PI = np.pi


def fit(input_data: list | np.ndarray, x: list | np.ndarray, fitmodel: Model, 
        fixed_parameters: list[str] = None, **fitting_kwargs):
    """
    Fit data based on a given mathematical model.
    User can choose Model in this file or built-in lmfit Model.
    Return to a ModelResult object.
    Fitting result can be accessed by result.best_values (dict).
    The corresponding fit-datapoint can be accessed by result.best_fit (ndarray).
    We allow keyword arguments here for 2D/3D fit and possible change.
    """
    input_data = np.array(input_data)
    x = np.array(x)
    fixed_parameters = make_it_list(fixed_parameters)

    fitmodel = fitmodel()
    params = fitmodel.guess(input_data, x, **fitting_kwargs)
    for param in fixed_parameters:
        params[param].vary = False
        
    return fitmodel.fit(input_data, params=params, x=x, **fitting_kwargs)


def exp_sin_func(x, tau, freq, phase, A, C):
    return C + A * exp(-x/tau) * sin(2*PI*freq*x + phase)

def exp_func(x, tau, A, C):
    return C + A * exp(-x/tau)

def quad_func(x, x0, A, C):
    return C + A * (x - x0)**2

def sin_func(x, freq, phase, A, C):
    return C + A * sin(2*PI*freq*x + phase)

def gaussian1d_func(x, mean, std, A, C):
    return C + A * exp( -0.5 * ((x-mean) / std)**2 )

def exp_func2(x, r, A, C):
    return C + A * (r ** x)


def chevron_func(x, y, omega_0, freq_offset, phase, A, C):
    """
    Ref: Gerry and Knight(2005) Eq.(4.80)
    Here x is the time/pulse_length, y is the detuning.
    """
    t, detuning = np.meshgrid(x, y)  # Return 2D data when t and detuning are both array.
    omega_R = np.sqrt( omega_0**2 + (detuning - freq_offset)**2 )
    return C + A * ( omega_0 / omega_R * sin(2*PI * omega_R * t / 2 + phase) )**2


def spectroscopy_func(x, t, omega_0, freq_offset, A, C):
    """
    Ref: Gerry and Knight(2005) Eq.(4.80)
    Here x is the detuning.
    """
    omega_R = np.sqrt( omega_0**2 + (x - freq_offset)**2 )
    return C + A * ( omega_0 / omega_R * sin(2*PI * omega_R * t / 2) ) ** 2


def resonator_hanger_transmission_func(x, f0, Q, Qc, theta, A, phi, ED, PCC):
    """
    Calculate S21 as a function of frequency.
    Here x is frequency, A and phi are global amplitude/phase factor.
    A and phi will effectively change the off-resonant amplitude/phase.
    ED is the electrical delay in unit of [s].
    PCC is Power Compression Coefficient to compensate natural power drop as mod_freq.
    PCC is in unit of [Amplitude/Hz], which is [mV/Hz] on Qblox.
    Reference:
    https://iopscience.iop.org/article/10.1088/2058-9565/ac070e
    """
    S21 =  1 - Q * exp(1j * theta) / Qc / (1 +  2j * Q * (x-f0) / f0)
    return A * exp(1j * (2*PI * (x-f0) * ED + phi)) * (1 + PCC * (x-f0)) * S21


def double_exp_sin_func(x, C_0, C_1, tau1, tau2R, A_0, A_1, freq_0, freq_1, phase_0, phase_1):
    """ 
    Ref: https://doi.org/10.1103/PhysRevLett.114.010501, Eq.(7) in Supplment Material.
    The function here is more general and can give better fit with T1 and separated amp.
    """
    return (C_0 + C_1 * exp(-x / tau1) 
            + exp(-x / tau2R) * (A_0 * sin(2*PI * freq_0 * x + phase_0) 
                                + A_1 * sin(2*PI * freq_1 * x + phase_1)))


def triple_exp_sin_func(x, C_0, C_1, tau1, tau2R, A_0, A_1, A_2, 
                        freq_0, freq_1, freq_2, phase_0, phase_1, phase_2):
    return (C_0 + C_1 * exp(-x / tau1) 
            + exp(-x / tau2R) * (A_0 * sin(2*PI * freq_0 * x + phase_0) 
                                + A_1 * sin(2*PI * freq_1 * x + phase_1)
                                + A_2 * sin(2*PI * freq_2 * x + phase_2)))


class ExpSinModel(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(func=exp_sin_func, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        sin_model = SineModel()
        sin_params = sin_model.guess(data, x=x)
        
        self.set_param_hint('freq', value=sin_params['frequency'].value/2/PI, min=0)
        self.set_param_hint('A', value=sin_params['amplitude'].value)
        self.set_param_hint('C', value=np.mean(data))
        self.set_param_hint('phase', value=sin_params['shift'].value)
        self.set_param_hint('tau', value=(x[-1] - x[0])/2, min=0)
        return self.make_params()      


class ExpModel(Model):
    """ Please do not use the built-in ExponentialModel of lmfit, because its guess is terrible.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=exp_func, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        self.set_param_hint('A', value=np.max(data) - np.min(data))
        self.set_param_hint('C', value=np.mean(data))
        self.set_param_hint('tau', value=(x[-1] - x[0])/2, min=0)
        return self.make_params()   
    
    
class QuadModel(Model):
    """ x0 is much more convenient than expression in  'a,b,c' form.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=quad_func, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        quad_model = QuadraticModel()
        quad_params = quad_model.guess(data, x=x)
        
        x0 = -quad_params['b']/quad_params['a']/2
        A = quad_params['a']
        C = quad_params['c'] - A * x0**2
        
        self.set_param_hint('x0', value=x0)
        self.set_param_hint('A', value=A)
        self.set_param_hint('C', value=C)
        return self.make_params() 
    
    
class SinModel(Model):
    """ Consider the 2pi problem and add offset, comparing to SineModel.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=sin_func, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        sin_model = SineModel()
        sin_params = sin_model.guess(data, x=x)
        
        self.set_param_hint('freq', value=sin_params['frequency'].value/2/PI, min=0)
        self.set_param_hint('A', value=sin_params['amplitude'].value)
        self.set_param_hint('C', value=np.mean(data))
        self.set_param_hint('phase', value=sin_params['shift'].value)
        return self.make_params()  
    
    
class ExpModel2(Model):
    """ A different exponential model where we fit base rather than life time.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=exp_func2, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        self.set_param_hint('A', value=np.max(data) - np.min(data))
        self.set_param_hint('C', value=np.min(data), min=0)
        self.set_param_hint('r', value=0.99, min=0)
        return self.make_params()
    

class ChevronModel(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(func=chevron_func, independent_vars=['x', 'y'], *args, **kwargs)
    
    def guess(self, data: np.ndarray, x: np.ndarray, **fitting_kwargs):
        self.set_param_hint('omega_0', value=2/(np.max(x)-np.min(x)))  # Assume two periods.
        self.set_param_hint('freq_offset', value=0)
        self.set_param_hint('phase', value=0)
        self.set_param_hint('A', value=np.max(data)-np.min(data))
        self.set_param_hint('C', value=np.mean(data))
        return self.make_params()
    

class SpectroscopyModel(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(func=spectroscopy_func, *args, **kwargs)
    
    def guess(self, data: np.ndarray, x: np.ndarray, t: float, **fitting_kwargs):
        # We make 't' as a fixed parameter since we don't need it.
        # The initial guess for freq_offset is crucial. It highly decides the fitting result.
        # We need to consider the peak can be both upward and downward.
        peak_idx = np.argmax(np.abs( data - np.mean(data) ))

        self.set_param_hint('t', value=t, vary=False)
        self.set_param_hint('omega_0', value=0.5 / t, min=0)
        self.set_param_hint('freq_offset', value=x[peak_idx], min=x[0], max=x[-1])
        self.set_param_hint('A', value=0, max=2*(np.max(data)-np.min(data)), min=-2*(np.max(data)-np.min(data)))
        self.set_param_hint('C', value=np.mean(data))
        return self.make_params()


class ResonatorHangerTransmissionModel(Model):
    def __init__(self, *args, **kwargs):
        super().__init__(func=resonator_hanger_transmission_func, *args, **kwargs)

    def guess(self, data: np.ndarray, x: np.ndarray, **fitting_kwargs):
        amplitude = np.abs(data)
        phase = np.angle(data)

        A_guess = np.max(amplitude)
        f0_guess = x[np.argmin(amplitude)]
        amplitude /= A_guess  # Normalize amplitude for later guess.

        # Guess all parameters except Qs.
        self.set_param_hint('f0', value=f0_guess, min=0)
        self.set_param_hint('theta', value=0, min=-PI, max=PI)
        self.set_param_hint('A', value=A_guess, min=0)
        self.set_param_hint('phi', value=np.mean(phase), min=-PI, max=PI)
        self.set_param_hint('ED', value=(phase[-1] - phase[0]) / (x[-1] - x[0]))
        self.set_param_hint('PCC', value=(amplitude[-1] - amplitude[0]) / (x[-1] - x[0]))

        # Guess Q.
        half_max = (np.max(amplitude) + np.min(amplitude)) / 2
        full_width = 2 * np.abs(f0_guess - x[np.argmin(np.abs(amplitude - half_max))])
        self.set_param_hint('Q', value=f0_guess/full_width, min=0)

        # Guess Qc.
        # Ref: https://lmfit.github.io/lmfit-py/examples/example_complex_resonator_model.html
        self.set_param_hint('Qc', value=f0_guess/full_width/(1-np.min(amplitude)), min=0)

        return self.make_params()  
    

class DoubleExpSinModel(Model):
    """ Designed to fit Ramsey experiment with two possible frequency adding together.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=double_exp_sin_func, *args, **kwargs)
        
    def guess(self, data: np.ndarray, x: np.ndarray, **fitting_kwargs):
        """
        fitting_kwargs take:
            height: float. Specify the minimal height to be identified as peak.
            fixed_params: dict. The parameters to fix during fitting.
                Example: {'phase_0': -PI/2, 'phase_1': -PI/2}
        """
        # Substract DC signal.
        data_mean = np.mean(data)
        data = data - data_mean  # Do not use -= or += here! It will change array in place.
        rfft_freqs = np.fft.rfftfreq(len(x), abs(x[-1] - x[0]) / (len(x) - 1))
        freq_step = rfft_freqs[1] - rfft_freqs[0]

        # DFT
        rfft = abs(np.fft.rfft(data))

        # Find peaks
        height = fitting_kwargs['height'] if 'height' in fitting_kwargs else 10 * max(rfft) / len(rfft)
        peaks, _ = find_peaks(rfft, height)
        if len(peaks) != 2: raise ValueError(f'fitting: Cannot fit {len(peaks)} peaks')

        # Set initial guess of parameters.
        f_0 = rfft_freqs[peaks[0]]
        f_1 = rfft_freqs[peaks[1]]
        A_0 = rfft[peaks[0]] / len(rfft)
        A_1 = rfft[peaks[1]] / len(rfft)
        tau = x[-1] - x[0]

        self.set_param_hint('freq_0', value=f_0, min=f_0-freq_step, max=f_0+freq_step)
        self.set_param_hint('freq_1', value=f_1, min=f_1-freq_step, max=f_1+freq_step)
        self.set_param_hint('A_0', value=A_0, min=A_0/10)
        self.set_param_hint('A_1', value=A_1, min=A_1/10)
        self.set_param_hint('C_0', value=data_mean/2)
        self.set_param_hint('C_1', value=data_mean/2)
        self.set_param_hint('tau1', value=tau, min=tau/10)
        self.set_param_hint('tau2R', value=tau, min=tau/10)
        self.set_param_hint('phase_0', value=0, min=-PI, max=PI)
        self.set_param_hint('phase_1', value=0, min=-PI, max=PI)

        if 'fixed_params' in fitting_kwargs:
            for param, value in fitting_kwargs['fixed_params'].items():
                self.set_param_hint(param, value=value, vary=False)

        return self.make_params() 
    

class TripleExpSinModel(Model):
    """ Designed to fit Ramsey experiment with three possible frequency adding together.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(func=triple_exp_sin_func, *args, **kwargs)
        
    def guess(self, data, x, **fitting_kwargs):
        # Substract DC signal.
        data_mean = np.mean(data)
        data = data - data_mean  # Do not use -= or += here! It will change array in place.
        rfft_freqs = np.fft.rfftfreq(len(x), abs(x[-1] - x[0]) / (len(x) - 1))
        freq_step = rfft_freqs[1] - rfft_freqs[0]

        # DFT
        rfft = abs(np.fft.rfft(data))

        # Find peaks
        height = fitting_kwargs['height'] if 'height' in fitting_kwargs else 15 * max(rfft) / len(rfft)
        peaks, _ = find_peaks(rfft, height)
        if len(peaks) != 3: raise ValueError(f'fitting: Cannot fit {len(peaks)} peaks')

        # Set initial guess of parameters.
        f_0 = rfft_freqs[peaks[0]]
        f_1 = rfft_freqs[peaks[1]]
        f_2 = rfft_freqs[peaks[2]]
        A_0 = rfft[peaks[0]] / len(rfft)
        A_1 = rfft[peaks[1]] / len(rfft)
        A_2 = rfft[peaks[2]] / len(rfft)
        tau = x[-1] - x[0]

        self.set_param_hint('freq_0', value=f_0, min=f_0-freq_step, max=f_0+freq_step)
        self.set_param_hint('freq_1', value=f_1, min=f_1-freq_step, max=f_1+freq_step)
        self.set_param_hint('freq_2', value=f_2, min=f_2-freq_step, max=f_2+freq_step)
        self.set_param_hint('A_0', value=A_0, min=A_0/10)
        self.set_param_hint('A_1', value=A_1, min=A_1/10)
        self.set_param_hint('A_2', value=A_2, min=A_2/10)
        self.set_param_hint('C_0', value=data_mean/2)
        self.set_param_hint('C_1', value=data_mean/2)
        self.set_param_hint('tau1', value=tau, min=tau/10)
        self.set_param_hint('tau2R', value=tau, min=tau/10)
        self.set_param_hint('phase_0', value=0, min=-PI, max=PI)
        self.set_param_hint('phase_1', value=0, min=-PI, max=PI)
        self.set_param_hint('phase_2', value=0, min=-PI, max=PI)

        if 'fixed_params' in fitting_kwargs:
            for param, value in fitting_kwargs['fixed_params'].items():
                self.set_param_hint(param, value=value, vary=False)

        return self.make_params() 