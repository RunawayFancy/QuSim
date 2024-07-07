import socket
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
import numpy as np

OPTIONS = {

}

class qsviewer:
    def __init__(self, name, port, host = 'localhost', *options):
        self.name = name
        self.port = port
        self.host = host
        if not options:
            self.opt = None
        else:
            self.opt = options
    
    def start_server(self):
        server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server_socket.bind((self.host, self.port))
        server_socket.listen(1)
        print(f"Listening on {self.host}:{self.port}...")

    def exec_app(self):
        app = pg.mkQApp(self.name)
        win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
        win.resize(1000,600)
        win.setWindowTitle('pyqtgraph example: Plotting')
        pg.setConfigOptions(antialias=True)
        pg.exec()

    def init_plot(self):
        global Plst

    def update(self):
        global curve, data, ptr, p6
        curve.setData(data[ptr%10])
        if ptr == 0:
            p6.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
        ptr += 1