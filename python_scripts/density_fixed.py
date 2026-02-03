import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
import matplotlib.animation as animation

# Check arguments
if len(sys.argv) < 4:
    print("Usage: python3 script.py <path> <x|y><index> <species1> [<species2> ...] [timestep]")
    sys.exit(1)

# Parse required arguments
file_name = 'result.h5'
path = sys.argv[1]
slice_arg = sys.argv[2].lower()

if not (slice_arg.startswith('x') or slice_arg.startswith('y')) or not slice_arg[1:].isdigit():
    print("Error: slice argument must be in format x<index> or y<index> (e.g., x10 or y50)")
    sys.exit(1)

slice_dir = slice_arg[0]
slice_index = int(slice_arg[1:])

# Parse species names and optional timestep
species_list = []
timestep_arg = None
for arg in sys.argv[3:]:
    if arg.isdigit():
        timestep_arg = int(arg)
    else:
        species_list.append(arg)

# Open HDF5 file
file_path = os.path.join(path, file_name)
f = h5py.File(file_path, 'r')
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DATA_TS_PHASE = int(NUM_TS / write_interval) + 1

# Set up plots
fig, axes = plt.subplots(2, 1, figsize=(12, 8), dpi=100, sharex=True)
ax_den, ax_pot = axes
colors = ['blue', 'green', 'orange', 'purple', 'magenta', 'cyan', 'brown']

def plot_frame(i):
    ax_den.clear()
    ax_pot.clear()

    for idx, species in enumerate(species_list):
        try:
            den2d = f[f"fielddata/den_{species}/{i}"][:]
        except KeyError:
            print(f"Warning: species '{species}' not found at timestep {i}. Skipping.")
            continue

        if slice_dir == "x":
            data_1d = den2d[slice_index, :]
            axis_label = f"x = {slice_index}"
        else:
            data_1d = den2d[:, slice_index]
            axis_label = f"y = {slice_index}"

        x_vals = np.arange(len(data_1d))
        ax_den.plot(x_vals, data_1d, label=species, color=colors[idx % len(colors)])

    # Potential
    pot2d = f[f"fielddata/pot/{i}"][:]
    pot_1d = pot2d[slice_index, :] if slice_dir == "x" else pot2d[:, slice_index]
    ax_pot.plot(x_vals, pot_1d, label='Potential', color='red')

    ax_den.set_ylabel(f"Density ({axis_label})")
    ax_pot.set_ylabel("Potential")
    ax_pot.set_xlabel("Grid Index")

    ax_den.set_title(f"Densities at timestep {i}")
    ax_pot.set_title(f"Potential at timestep {i}")

    ax_den.grid(True)
    ax_pot.grid(True)
    ax_den.legend()
    ax_pot.legend()

def on_key(event):
    if event.key == 'enter':
        on_key.frame = min(on_key.frame + 1, DATA_TS_PHASE - 1)
        plot_frame(on_key.frame * write_interval)
        plt.draw()
    elif event.key == 'backspace':
        on_key.frame = max(on_key.frame - 1, 0)
        plot_frame(on_key.frame * write_interval)
        plt.draw()
    elif event.key == 'tab':
        on_key.ani = animation.FuncAnimation(fig, lambda j: plot_frame(j * write_interval),
                                             frames=DATA_TS_PHASE, blit=False, interval=100, repeat=False)
        plt.draw()

on_key.frame = 0
on_key.ani = None

if timestep_arg is not None:
    if timestep_arg > NUM_TS:
        print(f"Error: timestep {timestep_arg} exceeds maximum {NUM_TS}")
        sys.exit(1)
    plot_frame(timestep_arg)
    plt.tight_layout()
    plt.show()
else:
    fig.canvas.mpl_connect('key_press_event', on_key)
    plot_frame(on_key.frame * write_interval)
    plt.tight_layout()
    plt.show()
