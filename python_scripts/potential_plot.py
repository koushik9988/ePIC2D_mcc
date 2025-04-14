import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py
import matplotlib.animation as animation

# hdf5 file name and path
file_name = 'result.h5'

path = sys.argv[1]

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DATA_TS_PHASE = int(NUM_TS / write_interval) + 1



#fig,ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10, 8), dpi=100)

def animate(j):
    i = j*write_interval*10
    #i = NUM_TS
    ax.clear()
    #pot = f[f"fielddata/pot/{i}"][:]
    pot = f[f"fielddata/den_electron/{i}"][:]
    #pot = f[f"fielddata/den_ion/{i}"][:]
    #pot = f[f"fielddata/efieldy/{i}"][:]
    pot = np.transpose(pot)
    cax = ax.imshow(pot, origin='lower', interpolation='bilinear', cmap='coolwarm')
    ax.set_title(f'Ion density at Timestep {i}')
    #fig.colorbar(cax, ax=ax, orientation='vertical')
    
    ax.set_aspect('equal')

    return cax

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
