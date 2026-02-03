import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

file_name = 'result.h5'
path = sys.argv[1]
output_dir = './plots'
path_fig = pjoin(path, output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Open HDF5 file
f = h5py.File(pjoin(path, file_name), 'r')

# --- Read number of species from species metadata section ---
spno = f["/metadata"].attrs.get("spno", None)
if spno is None:
    raise ValueError("Could not find 'spno' attribute in /metadata.")

print(f"Number of species : {spno}")

# --- Read species order from species metadata ---
species_order = f["/metadata_species"].attrs.get("species_order", None)
if species_order is None:
    raise ValueError("Could not find 'species_order' attribute in /metadata_species.")

print(f"Species order: {species_order}")

# Read kinetic energy matrix
data = f["time_var/kinetic_energy"]
ts = data[:, 0]  # First column is time

# Total number of plots: kinetic plots + potential + total
total_plots = spno + 2

# Calculate grid size (square-like)
cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

total_ke = np.zeros_like(ts)

# Plot kinetic energy for each species based on species_order
for i, species_name in enumerate(species_order):
    # Find index of the species in the file's metadata
    species_group = f["/metadata_species"].get(species_name, None)
    if species_group is None:
        raise ValueError(f"Species '{species_name}' not found in '/metadata_species'.")

    # Extract kinetic energy for the species
    start_idx = 1 + i * 3  # Each species has 3 KE components: x, y, z
    ke_x = data[:, start_idx]
    ke_y = data[:, start_idx + 1]
    ke_z = data[:, start_idx + 2]

    total_ke += ke_x + ke_y + ke_z  # Sum kinetic energies

    axs[i].plot(ts, ke_x, label=f"${species_name}~KE_x$")
    axs[i].plot(ts, ke_y, label=f"${species_name}~KE_y$")
    axs[i].plot(ts, ke_z, label=f"${species_name}~KE_z$")
    #axs[i].plot(ts, ke_x + ke_y + ke_z, label=f"${species_name}~KE_T$")
    axs[i].set_ylabel(f"{species_name} KE")
    axs[i].legend(loc='upper right', framealpha=0.5)

# Potential energy plot
pex = data[:, 3 * spno + 1]
pey = data[:, 3 * spno + 2]
potential_energy = pex + pey

axs[spno].plot(ts, pex, label="$\int E_x^2 dA$")
axs[spno].plot(ts, pey, label="$\int E_y^2 dA$")
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno].set_ylabel("Potential Energy")
axs[spno].legend(loc='upper right', framealpha=0.5)

# Total energy plot
total_energy = total_ke + potential_energy
axs[spno + 1].plot(ts, total_ke, label="KE", color='red')
axs[spno + 1].plot(ts, potential_energy, label="PE", color='green')
axs[spno + 1].plot(ts, total_energy, label="Total Energy", color='black')
#axs[spno + 1].set_ylim(np.min(total_energy) - 5, np.max(total_energy) + 5)
axs[spno + 1].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_ylabel("Total Energy")
axs[spno + 1].legend(loc='upper right', framealpha=0.5)

# Hide unused subplots if any
for i in range(total_plots, len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.show()
