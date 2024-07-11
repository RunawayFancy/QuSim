#!usr/bin/env python
# -*- coding: utf-8 -*-

from qusim.PulseGen.pulse_config import PulseConfig
from qusim.PulseGen.simulation_option import SimulationOption

import copy
import numpy as np; from numpy import pi as PI
import pyqtgraph as pg
from collections import defaultdict as ddict


def __pyqt_plot_test(pulse: PulseConfig, sim_opts: SimulationOption) -> None:
    '''
    do not use this in practice
    this is a tiny demo to show how pyqt works
    '''
    t_steps = sim_opts.tlist

    # waveform: 1d np array, e.g. (2400,)
    _waveform = pulse.get_pulse(sim_opts)
    _hint = (
        pulse.qindex,
        pulse.pulse_type
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


def plot_pulse_sequence(pseq: list[PulseConfig], sim_opts: SimulationOption) -> None:
    """
    Plot a sequence of pulses in a single QT window
    """

    t_steps = sim_opts.tlist
    pulse_map = ddict(list)

    for pulse in pseq:
        _waveform = pulse.get_pulse(sim_opts)/2/PI

        # _hint = (
        #     pulse.qindex,
        #     pulse.pulse_type,
        # )
        _hint = f"{pulse.pulse_type}{pulse.qindex}"
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
    from qusim.PulseGen.pulse_config import __TEST_PULSE__
    pulse = __TEST_PULSE__
    t_total = pulse.t_width + pulse.t_plateau + 20
    simopt = SimulationOption(
        simulation_time=t_total,
        simulation_step=int(100*t_total),
        initial_state=0
    )

    pseq = [copy.deepcopy(pulse) for _ in range(10)]
    ch_name = ['XY', 'Z', 'INT']
    for ii, pp in enumerate(pseq):
        pseq[ii].qindex = ii % 5
        pseq[ii].pulse_type = ch_name[ii % 3]
        pseq[ii].t_delay = ii*2

    plot_pulse_sequence(pseq, simopt)
