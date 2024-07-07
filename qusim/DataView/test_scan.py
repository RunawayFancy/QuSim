#%%
from qtplot import run_plotting_app
from multiprocessing import Process, Queue
import numpy as np
import time

def simulate_scan(data_queue):
    for i in range(50):
        data = np.random.normal(loc=0.0, scale=1.0, size=100)  # Generate some data
        data_queue.put(data)  # Push data to the queue
        time.sleep(0.1)  # Simulate time delay between data generations
#%%
data_queue = Queue()
plot_process = Process(target=run_plotting_app, args=(data_queue,))
scan_process = Process(target=simulate_scan, args=(data_queue,))

plot_process.start()
scan_process.start()

plot_process.join()
scan_process.join()

#%%

import pyqtgraph.examples
pyqtgraph.examples.run()

#%%

import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore

# Initialize the application or get the existing one
app = QtWidgets.QApplication.instance()
if app is None:
    app = QtWidgets.QApplication([])

# Use a global variable to keep track of the window
global win
try:
    # Try to update the existing window if it's already been created
    win.setWindowTitle('pyqtgraph example: Plotting - Updated')
except NameError:
    # If 'win' does not exist, create a new window
    win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
    win.resize(1000,600)
    win.setWindowTitle('pyqtgraph example: Plotting')

    pg.setConfigOptions(antialias=True)

    # Create plots as before
    p1 = win.addPlot(title="Basic array plotting", y=np.random.normal(size=100))
    p2 = win.addPlot(title="Multiple curves")
    p2.plot(np.random.normal(size=100), pen=(255,0,0), name="Red curve")
    p2.plot(np.random.normal(size=110)+5, pen=(0,255,0), name="Green curve")
    p2.plot(np.random.normal(size=120)+10, pen=(0,0,255), name="Blue curve")

    # Define the update function inside try-except block
    def update():
        global curve, data, ptr, p6
        curve.setData(data[ptr%10])
        if ptr == 0:
            p6.enableAutoRange('xy', False)  # Stop auto-scaling after the first data set is plotted
        ptr += 1

    # Adding an updating plot
    p6 = win.addPlot(title="Updating plot")
    curve = p6.plot(pen='y')
    data = np.random.normal(size=(10,1000))
    ptr = 0

    timer = QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(50)

if __name__ == '__main__':
    QtWidgets.QApplication.instance().exec_()

