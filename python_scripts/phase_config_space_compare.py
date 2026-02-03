import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys
import matplotlib.animation as animation

# 🔧 Accept multiple species separated by comma
if len(sys.argv) != 4:
    print("Usage: python3 script.py <path> <particle_type(s)> <plot_type>")
    print("plot_type options: x-vx, y-vy, x-y, vx-vy")
    sys.exit(1)

file_name = 'result.h5'
path = sys.argv[1]
particle_types = sys.argv[2].split(',')  # 🔧 Can be multiple, e.g. electron,ion
plot_type = sys.argv[3]

plot_path = './plots'
path_fig = pjoin(path, plot_path)
if not os.path.exists(path_fig):
    os.makedirs(path_fig)

f = h5py.File(pjoin(path, file_name), 'r')

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

fig, ax = plt.subplots()

colors = ['b', 'r', 'g', 'm', 'c']  # 🔧 Color cycle

def animate(i):
    j = i * write_interval_phase
    ax.clear()

    # 🔧 Plot each species
    for idx, species in enumerate(particle_types):
        color = colors[idx % len(colors)]

        data_phase_pos = f[f"particle_{species}/pos{j}"]
        x = data_phase_pos[:, 0]
        y = data_phase_pos[:, 1]

        data_phase_vel = f[f"particle_{species}/vel{j}"]
        vx = data_phase_vel[:, 0]
        vy = data_phase_vel[:, 1]

        if plot_type == 'x-vx':
            ax.scatter(x, vx, marker='.', color=color, alpha=0.8, s=11, label=f"{species}")
            ax.set_xlabel('$x$')
            ax.set_ylabel('$v_x$')
        elif plot_type == 'y-vy':
            ax.scatter(y, vy, marker='.', color=color, alpha=0.8, s=11, label=f"{species}")
            ax.set_xlabel('$y$')
            ax.set_ylabel('$v_y$')
        elif plot_type == 'x-y':
            ax.scatter(x, y, marker='.', color=color, alpha=0.8, s=11, label=f"{species}")
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
        elif plot_type == 'vx-vy':
            ax.scatter(vx, vy, marker='.', color=color, alpha=0.8, s=11, label=f"{species}")
            ax.set_xlabel('$v_x$')
            ax.set_ylabel('$v_y$')
        else:
            ax.text(0.5, 0.5, "Invalid plot type!", ha='center', va='center', transform=ax.transAxes)
            print(f"Invalid plot type: {plot_type}")

    ax.legend(loc='upper right', framealpha=0.5)
    return ax

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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None
fig.canvas.mpl_connect('key_press_event', on_key)

animate(on_key.frame)
plt.show()
