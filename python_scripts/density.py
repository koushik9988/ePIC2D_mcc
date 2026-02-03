import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 density.py <path> [species]")
    print("species defaults to the first entry in metadata_species")


if len(sys.argv) < 2:
    usage()
    sys.exit(1)

path = sys.argv[1]
file_name = "result.h5"

plot_dir = pjoin(path, "plots")
os.makedirs(plot_dir, exist_ok=True)

with h5py.File(pjoin(path, file_name), "r") as f:
    metadata = f["/metadata"]
    nx = int(metadata.attrs["Nx"])
    ny = int(metadata.attrs["Ny"])
    num_ts = int(metadata.attrs["NUM_TS"])
    write_int = int(metadata.attrs["write_int"])

    species_order = f["/metadata_species"].attrs.get("species_order", None)
    if species_order is None:
        print("Missing metadata_species/species_order")
        sys.exit(1)

    if len(sys.argv) >= 3:
        species = sys.argv[2]
    else:
        species = species_order[0]

    group_path = f"fielddata/den_{species}"
    if group_path not in f:
        print(f"Missing group: {group_path}")
        sys.exit(1)

    steps = list(range(0, num_ts + 1, write_int))
    if not steps:
        print("No timesteps found")
        sys.exit(1)

    sum_den = np.zeros((ny, nx))
    count = 0
    for ts in steps:
        ds_name = f"{group_path}/{ts}"
        if ds_name in f:
            sum_den += f[ds_name][:]
            count += 1

    if count == 0:
        print("No density datasets found")
        sys.exit(1)

    avg_den = sum_den / count

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        avg_den,
        origin="lower",
        extent=[0, nx, 0, ny],
        aspect="auto",
        cmap="viridis",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Time-averaged density: {species}")
    fig.colorbar(im, ax=ax, label="density")

    outname = f"avg_density_{species}.png"
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
