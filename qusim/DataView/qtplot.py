import sys
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg
from multiprocessing import Process, Queue
import numpy as np
import time

class RealTimePlot(QtWidgets.QWidget):
    def __init__(self, data_queue):
        super().__init__()
        self.data_queue = data_queue
        self.init_ui()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plot)
        self.timer.start(100)  # Update the plot every 100 ms

    def init_ui(self):
        self.graphWidget = pg.PlotWidget()
        self.plot = self.graphWidget.plot(pen=pg.mkPen(color=(255, 0, 0), width=2))

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.graphWidget)
        self.setLayout(layout)

        self.setGeometry(300, 300, 800, 600)
        self.setWindowTitle('Real-Time Plot')
        self.show()

    def update_plot(self):
        while not self.data_queue.empty():
            data = self.data_queue.get()
            self.plot.setData(data)  # Update the plot with the new data

def run_plotting_app(data_queue):
    app = QtWidgets.QApplication(sys.argv)
    window = RealTimePlot(data_queue)
    sys.exit(app.exec_())


