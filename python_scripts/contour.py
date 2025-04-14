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

# HDF5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

# Read HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DATA_TS_PHASE = int(NUM_TS / write_interval) + 1

# Set up figure and axis
fig, ax = plt.subplots(figsize=(10, 8), dpi=100)

def animate(j):
    i = j * write_interval * 10  # Timestep increment (adjustable)
    ax.clear()
    
    # Load data (e.g., electron density, ion density, potential, or efieldy)
    # Uncomment the desired field
    #pot = f[f"fielddata/pot/{i}"][:]         # Potential
    #pot = f[f"fielddata/den_electron/{i}"][:]  # Electron density
    #pot = f[f"fielddata/den_ion/{i}"][:]     # Ion density
    pot = f[f"fielddata/efieldy/{i}"][:]     # Y-component of electric field
    
    # Transpose data if needed (depends on your data orientation)
    pot = np.transpose(pot)
    
    # Create contour plot
    # Adjust levels (e.g., 20 contours) or specify custom levels if desired
    cax = ax.contourf(pot, levels=20, cmap='coolwarm')  # Filled contours
    contours = ax.contour(pot, levels=10, colors='black', linewidths=0.5)  # Line contours
    #ax.clabel(contours, inline=True, fontsize=8, fmt='%.2f')  # Label contours
    
    # Set title and aspect ratio
    #ax.set_title(f'Electron Density Contour at Timestep {i}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')
    
    # Add colorbar
    if not hasattr(animate, 'cbar'):  # Add colorbar only once
        animate.cbar = fig.colorbar(cax, ax=ax, orientation='vertical')
    else:
        animate.cbar.update_normal(cax)  # Update colorbar for new data
    
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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, 
                                             blit=False, interval=100, repeat=False)
        plt.draw()

# Initialize frame and animation reference
on_key.frame = 0
on_key.ani = None

# Connect key press event
fig.canvas.mpl_connect('key_press_event', on_key)

# Initial plot
animate(on_key.frame)
plt.show()

# Close HDF5 file (good practice)
f.close()