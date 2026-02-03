import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
import matplotlib.animation as animation

# Argument check
if len(sys.argv) < 3:
    print("Usage: python3 script.py <path> <field_type> [species_name]")
    print("field_type options: pot, efx, efy, density")
    sys.exit(1)

# Inputs
file_name = 'result.h5'
path = sys.argv[1]
field_type = sys.argv[2]

# Optional species name
species_name = sys.argv[3] if field_type == "density" and len(sys.argv) == 4 else None

# Open HDF5 file
f = h5py.File(os.path.join(path, file_name), 'r')
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DATA_TS_PHASE = int(NUM_TS / write_interval) + 1

# Optional grid mask overlay (1 = grid, 0 = empty)
grid_mask = metadata_group['grid_mask'][:] if 'grid_mask' in metadata_group else None

# Plot setup
fig, ax = plt.subplots(figsize=(10, 8), dpi=100)

def animate(j):
    i = j * write_interval * 10  
    ax.clear()

    if field_type == "pot":
        data = f[f"fielddata/pot/{i}"][:]
    elif field_type == "efx":
        data = f[f"fielddata/efieldx/{i}"][:]
    elif field_type == "efy":
        data = f[f"fielddata/efieldy/{i}"][:]
    elif field_type == "density":
        if species_name is None:
            print("Error: species_name argument required for density plot.")
            sys.exit(1)
        data = f[f"fielddata/den_{species_name}/{i}"][:]
    else:
        print(f"Unknown field_type: {field_type}")
        sys.exit(1)

    cax = ax.imshow(data, origin='lower', interpolation='bilinear', cmap='coolwarm', vmin=-10, vmax=10)
    if grid_mask is not None:
        ax.imshow(grid_mask, origin='lower', interpolation='nearest', cmap='Greys', alpha=0.35)
    ax.set_title(f'{field_type} at Timestep {i}')
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
        on_key.ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None

fig.canvas.mpl_connect('key_press_event', on_key)

animate(on_key.frame)
plt.show()
