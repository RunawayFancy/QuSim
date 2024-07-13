#!usr/bin/env python
# -*- coding: utf-8 -*-

from qusim.PulseGen.pulse_config import PulseConfig
from qusim.PulseGen.simulation_option import SimulationOption

import numpy as np
import pyqtgraph as pg
import socket
import pickle


class PGPlotServer():
    win: pg.GraphicsLayoutWidget
    s: socket.socket
    pseq: list[PulseConfig]
    
    def __init__(self, sim_opts: SimulationOption, host: str = '127.0.0.1', port: int = 63243, backlog: int = 10, refresh_msec: int = 500):
        self.win = pg.GraphicsLayoutWidget(show=True, title="Pulse Sequence")
        self.app = pg.mkQApp()
        self.pseq = []
        # x axis for the waveform, in np.ndarray
        self.sim_opts = sim_opts
        
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        self.s.bind((host, port))
        self.s.listen(backlog)
        
        timer = pg.QtCore.QTimer()
        timer.timeout.connect(self.auto_update)
        timer.start(refresh_msec)
        
        self.app.exec_()
    
    def __del__(self):
        self.s.close()
        
    def redraw(self):
        self.win.clear()
        p0 = None
        for _row, _pulse in enumerate(sorted(self.pseq)):
            _waveform = _pulse.get_pulse(self.sim_opts)/2/np.pi
            
            _plot = pg.PlotItem()
            self.win.addItem(_plot, row=_row, col=0)
        
            _plot.plot(x=self.sim_opts.tlist, y=_waveform)
            _plot.setLabel('left', text=f'{_pulse.pulse_type}{_pulse.qindex}')

            if p0:
                _plot.setXLink(p0)
                _plot.setYLink(p0)
            else:
                p0 = _plot
        
            if _row != len(self.pseq) - 1:
                # _plot.hideAxis('bottom')
                pass
            # show x axis lable, if bottom plot
            else:
                _plot.setLabel('bottom', text='Time (ns)')

    def auto_update(self):
        conn, _ = self.s.accept()
        with conn:
            while data := conn.recv(4096):
                new_pulse = pickle.loads(data)
                self.pseq.append(new_pulse)
                # print(new_pulse)

        # print('drawing') # trigger drawing
        self.redraw()
        

if __name__ == '__main__':
    simopt = SimulationOption(
        simu_time=24,
        simu_point=2400
    )
    pg_server = PGPlotServer(simopt)