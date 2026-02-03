import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 contour_xv.py <path> <species> <axis> [ts]")
    print("axis: x or y")


if len(sys.argv) < 4:
    usage()
    sys.exit(1)

path = sys.argv[1]
species = sys.argv[2]
axis = sys.argv[3].lower()
file_name = "result.h5"

if axis not in ("x", "y"):
    usage()
    sys.exit(1)

plot_dir = pjoin(path, "plots")
os.makedirs(plot_dir, exist_ok=True)

with h5py.File(pjoin(path, file_name), "r") as f:
    metadata = f["/metadata"]
    num_ts = int(metadata.attrs["NUM_TS"])
    write_int_phase = int(metadata.attrs["write_int_phase"])

    phase_group = f"particle_{species}"
    if phase_group not in f:
        print(f"Missing group: {phase_group}")
        sys.exit(1)

    steps = list(range(0, num_ts + 1, write_int_phase))
    if len(sys.argv) >= 5:
        ts = int(sys.argv[4])
        if f"pos{ts}" not in f[phase_group]:
            print(f"Timestep {ts} not found in {phase_group}")
            sys.exit(1)
    else:
        ts = steps[-1]

    pos = f[f"{phase_group}/pos{ts}"][:]
    vel = f[f"{phase_group}/vel{ts}"][:]

    axis_idx = 0 if axis == "x" else 1
    coord = pos[:, axis_idx]
    vcomp = vel[:, axis_idx]

    bins = 200
    hist, xedges, yedges = np.histogram2d(coord, vcomp, bins=bins)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        hist.T,
        origin="lower",
        aspect="auto",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="inferno",
    )
    ax.set_xlabel(axis)
    ax.set_ylabel(f"v{axis}")
    ax.set_title(f"{species}: {axis}-v{axis} at ts={ts}")
    fig.colorbar(im, ax=ax, label="count")

    outname = f"contour_{species}_{axis}v{axis}_ts{ts}.png"
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
