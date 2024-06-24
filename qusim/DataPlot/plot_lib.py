# -*- coding: utf-8 -*-
"""
@author: Jiheng Duan
"""
import qusim.PulseGen.pulse_waveform as pw
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def plot_pulse_sequence(pulse_sequence, simulation_option):
    t_list = np.linspace(0, simulation_option ["simulation_time"], 
     simulation_option ["simulation_step"])
    channel_dic = {}
    q_index_list = []
    for pulse in pulse_sequence:
        pulse_lib = pw.pulse_lib(pulse)
        pulse_amp = pulse["amplitude"]
        drive_pulse = pulse_lib.get_pulse(simulation_option)
        if pulse["q_index"] not in q_index_list:
            q_index_list.append(pulse["q_index"])
        # waveform_y = [drive_pulse(t, None) for t in t_list]  
        waveform_y = drive_pulse
        
        # Rescale the pulses by factor 1/1.2
        # for i in range(len(t_list)):
        #     if waveform_y[i] != None and pulse_amp != 0:
        #         waveform_y[i] /= np.abs(pulse_amp)

        channel_name = "{}{}".format(pulse["type"], pulse["q_index"])
        if channel_name in channel_dic:
            channel_dic[channel_name].append(waveform_y)
            # channel_dic[channel_name] = np.array(channel_dic[channel_name]) + np.array(waveform_y)
        else:
            channel_dic[channel_name] = [waveform_y]
    # print(type(channel_dic['XY0']))

    # Define a custom sorting function to sort the channel names as desired
    def custom_sort(channel_name):
        channel_type = channel_name[:-1]
        channel_index = channel_name[-1]
        if channel_index.isdigit():
            return int(channel_index), channel_type
        else:
            return float('inf'), channel_type

    # Sort the channel names based on the custom sorting function
    sorted_channels = sorted(channel_dic.keys(), key=custom_sort)

    # Convert the waveform data to NumPy arrays
    channel_dic = {k: np.array(v) for k, v in channel_dic.items()}

    # Create the figure and axes
    fig, ax1 = plt.subplots()
    fig.canvas.manager.set_window_title('pulse_sequence')
    # Set the vertical spacing between channel lines
    vertical_spacing = 2

    # Iterate over the channels and plot them vertically separated
    for i, channel_name in enumerate(sorted_channels):
        try:
            iterator = iter(channel_dic[channel_name])
        except TypeError:
            pass
        else:
            index = 1
            for waveform in channel_dic[channel_name]:
            # waveform = channel_dic[channel_name]
                vertical_offset = i * vertical_spacing
                waveform_offset = waveform + vertical_offset
                for jj in range(len(waveform_offset)):
                    if waveform_offset[jj] == vertical_offset: waveform_offset[jj] = None
                ax1.plot(t_list, waveform_offset, label=channel_name + str(index))
                index += 1

    # Set the y-axis tick labels and limits for channel plot
    ax1.set_yticks(np.arange(len(channel_dic)) * vertical_spacing)
    ax1.set_yticklabels(sorted_channels)
    ax1.set_ylim(-vertical_spacing, len(channel_dic) * vertical_spacing)
    ax1.set_xlim(0, simulation_option["simulation_time"])
    # Add legend and labels for channel plot
    # ax1.legend(loc="upper right")
    ax1.set_xlabel("Time/ns")
    ax1.set_ylabel("Channel")
    
    # Create a twin axes for the pulse amplitude
    ax2 = ax1.twinx()
    ax2.set_ylim(-vertical_spacing, len(channel_dic) * vertical_spacing)
    
    # Add label for the pulse amplitude axis
    ax2.set_ylabel("Normalized Pulse Amplitude")

    # Show the plot
    plt.grid()
    plt.show()

    return 0


def plot_population_evolution(_system, result_list, simulation_option, interested_state, interested_state_label, initial_state_label):
    num_subplots = len(result_list)
    # Calculate the number of rows and columns for the subplot layout
    num_rows = int(num_subplots ** 0.5)
    num_cols = (num_subplots + num_rows - 1) // num_rows
    # Create subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(10, 8))
    # Flatten the axs array if it's 1D (for the case when there's only one subplot)
    time_list = np.linspace(0,simulation_option["simulation_time"],simulation_option["simulation_step"])
    if num_subplots in [2,3]:
        axs = np.array([axs])
    if num_subplots == 1:
        axs = np.array([[axs]])
    for index, result in enumerate(result_list):
        row_idx = index // num_cols
        col_idx = index % num_cols
        label_list = interested_state_label[index]
        data_list = _system.get_data_list(result, simulation_option, interested_state[index])
        for jndex in range(len(data_list)):
            axs[row_idx, col_idx].plot(time_list, data_list[jndex], label = label_list[jndex])
            axs[row_idx, col_idx].set_xlabel("Time/ns")
            axs[row_idx, col_idx].set_ylabel("Population")
            axs[row_idx, col_idx].legend()
        axs[row_idx, col_idx].set_title("Initial state: {}".format(initial_state_label[index]))
    # Adjust layout
    plt.tight_layout()
    # Show the plot
    plt.show()
    return 0

def plot_zz_sweep(x_list:list, y_list:list, zz_list:list, x_label:str, y_label:str):
    nrm1 = matplotlib.colors.LogNorm(1e-6, 1e-1)  
    fig1, ax1 = plt.subplots()
    p1 = ax1.pcolor(x_list, y_list, zz_list, cmap=matplotlib.cm.viridis_r,norm=nrm1)
    cb1 = fig1.colorbar(p1, ax=ax1,label='ZZ/GHz', extend='both')
    ax1.set_xlabel('${}$/GHz'.format(x_label))
    ax1.set_ylabel('${}$/GHz'.format(y_label))

    plt.show()
    return 0

def plot_Elevel_dynamics(w_scan_space, energy_level_list, num_to_plot, label:str, xrange = [], yrange = []):
    if isinstance(num_to_plot, int):
        plot_list = range(num_to_plot)
    if isinstance(num_to_plot, list):
        if len(num_to_plot) ==2:
            if num_to_plot[1] - num_to_plot[0] != 1:
                plot_list = range(num_to_plot[0], num_to_plot[1])
            else:
                plot_list = num_to_plot
        else:
            plot_list = num_to_plot
    for j in plot_list:
        plt.plot(w_scan_space, [v[j] for v in energy_level_list])
    if len(xrange) > 1: plt.xlim(xrange)
    if len(yrange) > 1: plt.ylim(yrange)
    plt.xlabel('${}$/GHz'.format(label))
    plt.ylabel('Energy (GHz)')
    plt.title("Energy level")
    plt.show()
    return 0