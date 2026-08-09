import socket
import errno
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Set up the server socket
host = 'localhost'
port = 1234
server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server_socket.bind((host, port))
server_socket.listen(1)
print(f"Listening on {host}:{port}...")

# Prepare the plot
fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro-')

def init():
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 10)
    return ln,

def update(frame):
    if frame is None:
        return ln,  # No data to update, just return the line object

    # Store the current view limits before clearing
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Clear the previous data
    ax.cla()  # Clears the axes

    # Set new titles, labels, etc. here if necessary
    ax.set_title("Real-time Plot of Scan Data")
    ax.set_xlabel("X-axis Label")
    ax.set_ylabel("Y-axis Label")

    # Plot new data
    ax.plot(xdata, ydata, 'ro-')  # You can customize the plot style

    # Reapply the old view limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Redraw the canvas
    fig.canvas.draw()
    return ln,


def data_gen():
    conn, addr = server_socket.accept()
    conn.setblocking(0)  # Set socket to non-blocking mode
    print(f"Connected by {addr}")
    with conn:
        while True:
            try:
                data = conn.recv(1024)
                if not data:
                    # No data received, yield None to keep the loop alive
                    yield None
                    continue

                decoded_data = data.decode()
                x, y = map(float, decoded_data.split(','))
                xdata.append(x)
                ydata.append(y)
                ax.set_xlim(min(xdata), max(xdata) + 1)
                ax.set_ylim(min(ydata) - 1, max(ydata) + 1)
                yield x, y
            except socket.error as e:
                # Error because no data available
                if e.errno != errno.EWOULDBLOCK:
                    print(f"Socket error: {e}")
                    break
                yield None  # No data available, yield None to maintain the plot update
                continue


ani = FuncAnimation(fig, update, frames=data_gen, init_func=init, blit=True)
plt.show()


