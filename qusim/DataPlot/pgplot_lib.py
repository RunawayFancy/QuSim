#!usr/bin/env python
# -*- coding: utf-8 -*-
#%%
import sys
sys.path.append("../..")

from qusim.PulseGen.pulse_config import PulseConfig
from qusim.PulseGen.simulation_option import SimulationOption

from matplotlib.backends.backend_qtagg import FigureCanvas

import copy
import numpy as np; from numpy import pi as PI
import pyqtgraph as pg
from collections import defaultdict as ddict


class QTWaveformPlotter():
    win: pg.GraphicsLayoutWidget
    sim_opts: SimulationOption
    # _key <=> _plot <=> _waveform is 1-1-1 mapping
    _curve_hashmap: dict[PulseConfig, pg.PlotDataItem]
    _plotlast: pg.PlotItem | None
    _plot0: pg.PlotItem | None
    
    def __init__(self, sim_opts: SimulationOption, title: str = 'Pulse Sequence'):
        assert sim_opts is not None
        
        # a layout contains several plots
        self.win = pg.GraphicsLayoutWidget(show=True, title=title)
        self.sim_opts = sim_opts
        self._curve_hashmap = {}
        self._plot0 = None
        self._plotlast = None
        self.app = pg.mkQApp()
        
        # timer = pg.QtCore.QTimer()
        # timer.timeout.connect(update)
        # timer.start(50)

    
    def draw(self, *pseq: PulseConfig) -> None:
        # always plot new pulses below old ones
        for _pulse in sorted(pseq):
            # form key for searching two hashmaps
            # add waveform to waveform hashmap
            _val = _pulse.get_pulse(self.sim_opts)/2/PI
            qindex = _pulse.qindex

            # based on _pulse, add/find the plot
            if _pulse in self._curve_hashmap:
                _curve = self._curve_hashmap[_pulse]
                _curve.setData(x=self.sim_opts.tlist, y=_val)
            else:
                _plot = pg.PlotItem()
                self.win.addItem(_plot, row=len(self._curve_hashmap), col=0)
                _curve = _plot.plot(x=self.sim_opts.tlist, y=_val)
        
                self._curve_hashmap[_pulse] = _curve

                # add label and sync X axis or Y axis
                _plot.setLabel('left', text=f'{_pulse.pulse_type}{qindex}')
                # _plot.hideAxis('bottom')
                if self._plot0 is None:
                    self._plot0 = _plot
                else:
                    _plot.setXLink(self._plot0)
                    _plot.setYLink(self._plot0)
        
        # self.win.show()
        # self.app.exec_()


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
    pulse_map: ddict[PulseConfig, list] = ddict(list)

    for pulse in pseq:
        _waveform = pulse.get_pulse(sim_opts) / 2 / PI
        pulse_map[pulse].append(_waveform)

    win = pg.GraphicsLayoutWidget(show=True, title="pulse_sequence")
    p0 = None

    # reverse=True becuase QT plot from top to bottom
    # matplotlib plot from bottom to top along y-pos direction
    # so reverse the order to keep figures the same
    for _row, _pulse in enumerate(sorted(pulse_map.keys(), reverse=True)):
        _val = pulse_map[_pulse]

        _plot = pg.PlotItem()
        win.addItem(_plot, row=_row, col=0)

        _plot.multiDataPlot(x=[t_steps] * len(_val), y=_val)
        _plot.setLabel('left', text=f'{_pulse.pulse_type}{_pulse.qindex}')

        # link x axis to the first plot, if there is first plot
        if p0:
            _plot.setXLink(p0)
            _plot.setYLink(p0)
        else:
            p0 = _plot

        # hide x axis, if not bottom plot
        if _row != len(pulse_map) - 1:
            _plot.hideAxis('bottom')
        # show x axis lable, if bottom plot
        else:
            _plot.setLabel('bottom', text='Time (ns)')

    pg.mkQApp().exec_()


if __name__ == '__main__':
    from qusim.PulseGen.pulse_config import __TEST_PULSE__
    
    pulse = __TEST_PULSE__
    t_total = pulse.t_width + pulse.t_plateau + 20
    simopt = SimulationOption(
        simu_time=t_total,
        simu_point=int(100*t_total)
    )

    pseq = [copy.deepcopy(pulse) for _ in range(10)]
    ch_name = ['XY', 'Z', 'INT']
    for ii, pp in enumerate(pseq):
        pseq[ii].qindex = ii % 5
        pseq[ii].pulse_type = ch_name[ii % 3]
        pseq[ii].t_delay = ii*2

    plotter = QTWaveformPlotter(simopt)

    FigureCanvas(plotter)
    for i, p in enumerate(pseq):
        print(i)
        plotter.draw(p)
    # plot_pulse_sequence(pseq, simopt)
