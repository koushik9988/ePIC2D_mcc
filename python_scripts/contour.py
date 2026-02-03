import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 contour.py <path> <quantity> [ts]")
    print("quantity: pot | efieldx | efieldy | den_<species>")


if len(sys.argv) < 3:
    usage()
    sys.exit(1)

path = sys.argv[1]
quantity = sys.argv[2]
file_name = "result.h5"

plot_dir = pjoin(path, "plots")
os.makedirs(plot_dir, exist_ok=True)

with h5py.File(pjoin(path, file_name), "r") as f:
    metadata = f["/metadata"]
    nx = int(metadata.attrs["Nx"])
    ny = int(metadata.attrs["Ny"])

    if quantity in ("pot", "efieldx", "efieldy"):
        group_path = f"fielddata/{quantity}"
    elif quantity.startswith("den_"):
        group_path = f"fielddata/{quantity}"
    else:
        print(f"Unknown quantity: {quantity}")
        usage()
        sys.exit(1)

    if group_path not in f:
        print(f"Missing group: {group_path}")
        sys.exit(1)

    group = f[group_path]
    steps = sorted([int(k) for k in group.keys()])
    if not steps:
        print(f"No datasets found in {group_path}")
        sys.exit(1)

    if len(sys.argv) >= 4:
        ts = int(sys.argv[3])
        if str(ts) not in group:
            print(f"Timestep {ts} not found in {group_path}")
            sys.exit(1)
    else:
        ts = steps[-1]

    data = group[str(ts)][:]

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        data,
        origin="lower",
        extent=[0, nx, 0, ny],
        aspect="auto",
        cmap="plasma",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{quantity} at ts={ts}")
    fig.colorbar(im, ax=ax, label=quantity)

    outname = f"contour_{quantity}_ts{ts}.png"
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
