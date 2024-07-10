#!usr/bin/env python
# -*- coding: utf-8 -*-

from qusim.PulseGen.pulse_waveform import PulseClass

import copy
import numpy as np
import pyqtgraph as pg
from collections import defaultdict as ddict


def __pyqt_plot_test(pulse_config, **sim_options) -> None:
    '''
    do not use this in practice
    this is a tiny demo to show how pyqt works
    '''
    sim_duration = sim_options['simulation_time']
    sim_step = sim_options['simulation_step']
    initial_state = sim_options['initial_state']

    t_steps = np.linspace(0, sim_duration, sim_step)

    # waveform: 1d np array, e.g. (2400,)
    _waveform = PulseClass(pulse_config).get_pulse(sim_options)
    _hint = (
        pulse_config['q_index'],
        pulse_config['type']
    )

    win = pg.GraphicsLayoutWidget(show=True, title="window title here")
    p1 = win.addPlot(x=t_steps, y=_waveform, row=0, col=0, title='Pulse 1')
    p2 = win.addPlot(x=t_steps, y=_waveform, row=1, col=0, title='Pulse 2')

    p2.setXLink(p1)
    # v = w.addViewBox(row=1, col=0, colspan=2)

    # p1.plot(_waveform, pen='r')
    # p2.plot(_waveform, pen='r')
    pg.mkQApp().exec_()
    return


def plot_pulse_sequence(pulse_sequence: list, **sim_options) -> None:
    """
    """
    assert 'simulation_time' in sim_options
    assert 'simulation_step' in sim_options
    assert 'initial_state' in sim_options

    sim_duration = sim_options['simulation_time']
    sim_step = sim_options['simulation_step']
    initial_state = sim_options['initial_state']

    t_steps = np.linspace(0, sim_duration, sim_step)

    pulse_map = ddict(list)

    for pulse_config in pulse_sequence:
        _waveform = PulseClass(pulse_config).get_pulse(sim_options)

        _hint = (
            pulse_config['q_index'],
            pulse_config['type'],
        )

        pulse_map[_hint].append(_waveform)

    win = pg.GraphicsLayoutWidget(show=True, title="window title here")
    p0 = None

    # reverse=True becuase QT plot from top to bottom
    # matplotlib plot from bottom to top along y-pos direction
    # so reverse the order to keep figures the same
    for _row, _k in enumerate(sorted(pulse_map.keys(), reverse=True)):
        _v = pulse_map[_k]

        _p = pg.PlotItem()
        win.addItem(_p, row=_row, col=0)

        _p.multiDataPlot(x=[t_steps] * len(_v), y=_v)
        _p.setLabel('left', text=f'{_k[1]} {_k[0]}')

        # link x axis to the first plot, if there is first plot
        if p0:
            _p.setXLink(p0)
            _p.setYLink(p0)
        else:
            p0 = _p

        # hide x axis, if not bottom plot
        if _row != len(pulse_map) - 1:
            _p.hideAxis('bottom')
        # show x axis lable, if bottom plot
        else:
            _p.setLabel('bottom', text='Time (ns)')

    pg.mkQApp().exec_()


if __name__ == '__main__':
    pulse = {
        'pulse_index': 0,
        'type': "INT",
        'pulse_shape': "cosine",
        't_delay': 0,
        't_width': 4,
        't_plateau': 0,
        'freq': 3,
        'phase': 0,
        'amplitude': 0.5,
        'q_index': 0
    }
    pulse2 = copy.deepcopy(pulse)
    pulse2['t_delay'] = 2
    t_total = pulse['t_width'] + pulse['t_plateau'] + 20
    simopt = {"simulation_time": t_total,
              "simulation_step": int(100*t_total), "initial_state": 0}

    pseq = [copy.deepcopy(pulse) for _ in range(10)]
    ch_name = ['XY', 'Z', 'INT']
    for ii, pp in enumerate(pseq):
        pseq[ii]['q_index'] = ii % 5
        pseq[ii]['type'] = ch_name[ii % 3]
        pseq[ii]['t_delay'] = ii*2

    plot_pulse_sequence(pseq + [pulse2], **simopt)

    # plib = PulseClass(pulse)
    # tlst = np.linspace(0, simopt["simulation_time"], simopt["simulation_step"])
    # plt.plot(tlst, plib.get_pulse(simopt))
