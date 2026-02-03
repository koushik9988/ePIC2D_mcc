import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin

# Argument check
if len(sys.argv) != 4:
    print("Usage: python3 plot_single_particle.py <path> <particle_type> <particle_index>")
    sys.exit(1)

# Inputs
path = sys.argv[1]
particle_type = sys.argv[2]
particle_index = int(sys.argv[3])
file_name = 'result.h5'

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')

# Read metadata
metadata = f['/metadata']
NUM_TS = metadata.attrs['NUM_TS']
write_interval_phase = metadata.attrs['write_int_phase']
DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Trajectory arrays
x_traj = []
timesteps = []

# Loop over time steps and record particle position
for i in range(DATA_TS_PHASE):
    j = i * write_interval_phase

    pos_key = f"particle_{particle_type}/pos{j}"
    if pos_key not in f:
        print(f"Missing dataset: {pos_key}")
        continue

    pos = f[pos_key]

    if particle_index >= pos.shape[0]:
        print(f"Error: particle index {particle_index} out of range at time {j}")
        sys.exit(1)

    x = pos[particle_index, 0]  # x-component of position
    x_traj.append(x)
    timesteps.append(j)

# Plot x vs t
plt.figure(figsize=(8, 5))
plt.plot(timesteps, x_traj, marker='o', markersize=4, linewidth=1.5, label=f'Particle {particle_index}')
plt.xlabel("Time step")
plt.ylabel("x position")
plt.title(f"Trajectory of Particle {particle_index} ({particle_type})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
