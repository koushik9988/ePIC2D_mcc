import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from os.path import join as pjoin
import sys
import h5py

# HDF5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DATA_TS_PHASE = int(NUM_TS / write_interval) + 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def animate(j):
    i = j * write_interval
    #i = NUM_TS
    ax.clear()
    pot = f[f"fielddata/pot/{i}"][:]
    #pot = f[f"fielddata/efieldx/{i}"][:]
    x = np.arange(pot.shape[1])
    y = np.arange(pot.shape[0])
    x, y = np.meshgrid(x, y)
    ax.plot_surface(x, y, pot, cmap='viridis')
    ax.set_title(f'Potential Field at Timestep {i}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Potential')
    return ax,

ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=10, repeat=False)

plt.show()
