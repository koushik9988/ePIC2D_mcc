import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# Argument
if len(sys.argv) != 3:
    print("Usage: python3 script.py <path> <particle_type>")
    sys.exit(1)

# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]
particle_type = sys.argv[2]

plot_path = './plots'

path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
#metadata_group = f['/metadata']

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')


# Read attributes
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1


# Debye length calculation
ni0 = 1e13
ne0 = ni0


# Plotting
#fig, (ax1, ax2) = plt.subplots(2, 1)
fig, ax1 = plt.subplots()

def animate(i):
    j = i * write_interval_phase

    # Phase space data
    data_phase_pos = f[f"particle_{particle_type}/pos{j}"]
    x = data_phase_pos[:, 0]
    y = data_phase_pos[:, 1]

    data_phase_vel = f[f"particle_{particle_type}/vel{j}"]
    vx = data_phase_vel[:, 0]
    vy = data_phase_vel[:, 1]

   
    ax1.clear()
    ax1.scatter(y, vy, marker='.', color='b', alpha=1.0, s=11, label=f"{particle_type} config space")
    ax1.set_xlabel('$x$')
    ax1.set_ylabel('$v$')
    ax1.legend(loc='upper right', framealpha=0.5)
    #ax1.set_aspect('equal', adjustable='box')  # Set equal aspect ratio
    """
    ax2.clear()
    ax2.scatter(y, vy, marker='.', color='b', alpha=1.0, s=11, label=f"{particle_type} velocity space")
    #ax2.set_aspect('equal', adjustable='box')  # Set equal aspect ratio
    """

    return ax1

def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS_PHASE - 1)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        animate(on_key.frame)
        plt.draw()
    elif event.key == 'tab':
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval= 100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None  # To keep reference to animation

fig.canvas.mpl_connect('key_press_event', on_key)

animate(on_key.frame)
plt.show()
