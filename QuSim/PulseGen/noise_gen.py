# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import numpy as np
import numpy.fft as fft

def noise_gen(tstart, tstop, tstep, config):
    """
    Generate noise to be added to a qubit control pulse.

    Parameters:
    - tstart: Start time of the simulation.
    - tstop: Stop time of the simulation.
    - tstep: Time step of the simulation.
    - config: Dictionary containing noise configuration. 
              It includes noise type, parameters, noise starting and stopping time.

    Returns:
    - noise: A numpy array containing the generated noise over time.
    """
    # Time array
    t = np.linspace(tstart, tstop, tstep)
    dt = t[1] - t[0]
    noise = np.zeros_like(t)
    noise_base = np.ones_like(t)

    switch = config.get('switch', 'off')    
    if isinstance(switch, str):
        if switch in ['on', 'off']:
            if switch == 'off': return noise
        else: raise(ValueError('Invalid switch: {switch}, str(on) and str(off) only.'))
    else: raise(ValueError('Invalid switch type: {switch}, str(on) and str(off) only.'))
    
    # Check for the type of noise and generate accordingly
    noise_type = config.get('type', '').lower()
    tau0 = config.get('tau0', tstop)
    noise_segments = segmentize(t, tau0, noise_base)
    
    if noise_type == 'gaussian':
        mean = config.get('mean', 0)
        std = config.get('std', 1)
        for i in range(len(noise_segments)):
            noise_segments[i] *= np.random.normal(mean, std)
        noise = np.concatenate(noise_segments)

    elif noise_type == 'white':
        """
        White noise is Gaussian noise with a flat spectral density
        """ 
        std = config.get('std', 1)
        for i in range(len(noise_segments)):
            noise_segments[i] *= np.random.normal(0, std)
        noise = np.concatenate(noise_segments)
        
    elif noise_type == 'rt':
        """
        Random telegraphic noise alternates between two values at random intervals
        """
        high = config.get('high', 1)
        low = config.get('low', -1)
        p_switch = config.get('p_switch', 0.01)  # Probability to switch value
        value = low
        for i in range(len(noise)):
            noise[i] = value
            if np.random.rand() < p_switch:
                value = high if value == low else low

    elif noise_type == 'jn': 
        """
        Johnson Nyquiest noise
        """
        kB = 1.38e-23  # Boltzmann constant in J/K
        temperature = config.get('temperature', 300)  # Room temperature default
        resistance = config.get('resistance', 1)  # Resistance in Ohms
        
        # Johnson-Nyquist noise voltage (RMS)
        V_noise_rms = np.sqrt(4 * kB * temperature * resistance)
        
        # Since Johnson noise is white, we model it as Gaussian noise with a standard deviation of V_noise_rms
        for i in range(len(noise_segments)):
            noise_segments[i] *= np.random.normal(0, V_noise_rms)
        noise = np.concatenate(noise_segments)

    elif noise_type == '1/f':
        alpha = config.get('alpha', 1.0)  # Exponent in the 1/f^alpha relationship
        scale = config.get('scale', 1.0)  # Scale to control the noise amplitude

        # Create white noise as a base
        for i in range(len(noise_segments)):
            noise_segments[i] *= np.random.normal(0, 1)
        white_noise = np.concatenate(noise_segments)

        # Perform FFT
        f_transform = fft.fft(white_noise)

        # Generate frequency axis and avoid division by zero
        freq = fft.fftfreq(t.size, d=tstep)
        freq[0] = 1  # To avoid division by zero at the zero frequency

        # Adjust amplitude by 1/f^alpha
        f_adjusted = f_transform / np.abs(freq) ** (alpha / 2.0)

        # Inverse FFT to get back to time domain
        pink_noise = np.real(fft.ifft(f_adjusted))

        # Normalize to standard normal distribution
        pink_noise = (pink_noise - np.mean(pink_noise)) / np.std(pink_noise)

        # Apply the scale parameter to control the noise amplitude
        noise = pink_noise * scale

    elif noise_type == 'bistable':
        # Bistable noise alternates between two values with non-Markovian switching characteristics
        high = config.get('high', 1)
        low = config.get('low', -1)
        rate = config.get('rate', 0.1)  # Average rate of switching
        state_duration = 0  # Track duration in the current state
        
        # Initial state
        value = low if np.random.rand() < 0.5 else high
        noise[0] = value
        
        for i in range(1, len(noise)):
            noise[i] = value
            # Increase state duration
            state_duration += dt
            
            # Switching probability can depend on the state duration to introduce memory effect
            # For simplicity, we use a fixed rate here
            if np.random.rand() < rate * dt:
                value = high if value == low else low
                state_duration = 0  # Reset state duration after switching
    # Other noise types can be added here following a similar pattern
    # Apply noise only within the specified start and stop times, if given
    noise_start = config.get('noise_start', tstart)
    noise_stop = config.get('noise_stop', tstop)
    noise[(t < noise_start) | (t > noise_stop)] = 0
    
    return noise

def segmentize(t, tau0, noise):
    segments = []
    start_idx = 0

    while start_idx < len(t):
        # Find the end index for the current segment
        end_idx = start_idx
        while end_idx < len(t) and (t[end_idx] - t[start_idx]) <= tau0:
            end_idx += 1
        
        # Add the segment to the list
        segments.append(noise[start_idx:end_idx])
        
        # Move to the next starting index
        start_idx = end_idx
    
    return segments

# Todo list:
# 1. Add multi-level rt noise
# 2. Hyperfine noise