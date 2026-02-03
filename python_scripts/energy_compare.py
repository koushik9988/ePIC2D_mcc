import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import sys
from os.path import join as pjoin

# --- Setup ---
if len(sys.argv) != 3:
    raise ValueError("Usage: python compare_energy.py <path1> <path2>")

paths = [sys.argv[1], sys.argv[2]]
labels = ["Run 1", "Run 2"]
styles = ['-', '--']  # Line styles for runs

file_name = 'result.h5'
output_dir = './plots'
path_fig = pjoin(paths[0], output_dir)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# --- Open both files ---
files = [h5py.File(pjoin(p, file_name), 'r') for p in paths]

# --- Metadata ---
spno = files[0]["/metadata"].attrs.get("spno", None)
species_order = files[0]["/metadata_species"].attrs.get("species_order", None)

if spno is None or species_order is None:
    raise ValueError("Missing 'spno' or 'species_order' attribute in metadata.")

print(f"Number of species : {spno}")
print(f"Species order     : {species_order}")

# --- Setup plotting grid ---
total_plots = spno + 2
cols = int(np.ceil(np.sqrt(total_plots)))
rows = int(np.ceil(total_plots / cols))

fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True)
axs = axs.flatten()

# --- Color map for KE components ---
ke_colors = {'x': 'blue', 'y': 'orange', 'z': 'purple'}
pe_colors = {'Ex2': 'green', 'Ey2': 'magenta'}
total_colors = {'KE': 'red', 'PE': 'green', 'Total': 'black'}

# --- Plot ---
for idx, f in enumerate(files):
    data = f["time_var/kinetic_energy"]
    ts = data[:, 0]
    total_ke = np.zeros_like(ts)

    for i, species_name in enumerate(species_order):
        start_idx = 1 + i * 3
        ke_x = data[:, start_idx]
        ke_y = data[:, start_idx + 1]
        ke_z = data[:, start_idx + 2]

        total_ke += ke_x + ke_y + ke_z

        axs[i].plot(ts, ke_x, linestyle=styles[idx], color=ke_colors['x'], label=f"{labels[idx]}: KE_x")
        axs[i].plot(ts, ke_y, linestyle=styles[idx], color=ke_colors['y'], label=f"{labels[idx]}: KE_y")
        axs[i].plot(ts, ke_z, linestyle=styles[idx], color=ke_colors['z'], label=f"{labels[idx]}: KE_z")
        axs[i].set_ylabel(f"{species_name} KE")
        axs[i].legend(loc='upper right', framealpha=0.5)

    # --- Potential Energy ---
    pex = data[:, 3 * spno + 1]
    pey = data[:, 3 * spno + 2]
    potential_energy = pex + pey

    axs[spno].plot(ts, pex, linestyle=styles[idx], color=pe_colors['Ex2'], label=f"{labels[idx]}: $\int E_x^2$")
    axs[spno].plot(ts, pey, linestyle=styles[idx], color=pe_colors['Ey2'], label=f"{labels[idx]}: $\int E_y^2$")
    axs[spno].set_ylabel("Potential Energy")
    axs[spno].legend(loc='upper right', framealpha=0.5)

    # --- Total Energy ---
    total_energy = total_ke + potential_energy
    #axs[spno + 1].plot(ts, total_ke, linestyle=styles[idx], color=total_colors['KE'], label=f"{labels[idx]}: KE")
    #axs[spno + 1].plot(ts, potential_energy, linestyle=styles[idx], color=total_colors['PE'], label=f"{labels[idx]}: PE")
    axs[spno + 1].plot(ts, total_energy, linestyle=styles[idx], color=total_colors['Total'], label=f"{labels[idx]}: Total")
    axs[spno + 1].set_ylim(np.min(total_energy) - 5, np.max(total_energy) + 5)
    axs[spno + 1].set_ylabel("Total Energy")
    axs[spno + 1].legend(loc='upper right', framealpha=0.5)

# --- Label x-axis for bottom plots ---
axs[spno].set_xlabel("$\omega_{pe}t$")
axs[spno + 1].set_xlabel("$\omega_{pe}t$")

# --- Clean extra subplots ---
for i in range(total_plots, len(axs)):
    fig.delaxes(axs[i])

plt.tight_layout()
plt.show()
