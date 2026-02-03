import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.constants import value as constants
from os.path import join as pjoin
import os
import sys

# Argument check
if len(sys.argv) != 4:
    print("Usage: python3 plot_phase_time.py <path> <particle_type> <plot_type>")
    print("plot_type options: x-t, y-t, vx-t, vy-t, vz-t")
    sys.exit(1)

# File and plot setup
file_name = 'result.h5'
path = sys.argv[1]
particle_type = sys.argv[2]
plot_type = sys.argv[3]
plot_path = './plots'
path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read file
f = h5py.File(pjoin(path, file_name), 'r')

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read metadata
metadata_group = f['/metadata']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval_phase = metadata_group.attrs['write_int_phase']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Select variable to plot
var_map = {
    'x-t': ('pos', 0, 'x', '$x$'),
    'y-t': ('pos', 1, 'y', '$y$'),
    'vx-t': ('vel', 0, 'v_x', '$v_x$'),
    'vy-t': ('vel', 1, 'v_y', '$v_y$'),
    'vz-t': ('vel', 2, 'v_z', '$v_z$'),
}

if plot_type not in var_map:
    print(f"Invalid plot type: {plot_type}")
    sys.exit(1)

group_type, col_idx, axis_label, long_label = var_map[plot_type]

# Define bins
n_bins_pos = 200
n_bins_t = DATA_TS_PHASE
value_bins = None

# First pass: determine bin range
all_values = []
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data = f[f"particle_{particle_type}/{group_type}{j}"][:, col_idx]
    all_values.append(data)
all_values = np.concatenate(all_values)
vmin, vmax = np.percentile(all_values, [1, 99])  # robust bounds
value_bins = np.linspace(vmin, vmax, n_bins_pos + 1)

# 2D histogram initialization
hist2d = np.zeros((n_bins_pos, n_bins_t), dtype=np.int32)

# Fill 2D histogram (value vs time)
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase
    data = f[f"particle_{particle_type}/{group_type}{j}"][:, col_idx]
    hist, _ = np.histogram(data, bins=value_bins)
    hist2d[:, i] = hist

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
extent = [0, NUM_TS, value_bins[0], value_bins[-1]]
aspect = 'auto'

im = ax.imshow(hist2d, extent=extent, origin='lower', interpolation='bilinear', aspect=aspect, cmap='inferno')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Particle count')

ax.set_xlabel('Time step')
ax.set_ylabel(long_label)
ax.set_title(f"{long_label} vs Time ({particle_type})")

# Save and show
outname = f"{plot_type}_{particle_type}.png"
plt.savefig(pjoin(path_fig, outname), dpi=300)
print(f"Saved plot to {pjoin(path_fig, outname)}")
plt.show()
