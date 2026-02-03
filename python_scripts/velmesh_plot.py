import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 velmesh_plot.py <path> <species> [ts]")


if len(sys.argv) < 3:
    usage()
    sys.exit(1)

path = sys.argv[1]
species = sys.argv[2]
file_name = "result.h5"

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

    if len(sys.argv) >= 4:
        ts = int(sys.argv[3])
        if f"vel{ts}" not in f[phase_group]:
            print(f"Timestep {ts} not found in {phase_group}")
            sys.exit(1)
    else:
        ts = num_ts - (num_ts % write_int_phase)

    vel = f[f"{phase_group}/vel{ts}"][:]
    vx = vel[:, 0]
    vy = vel[:, 1]

    bins = 200
    hist, xedges, yedges = np.histogram2d(vx, vy, bins=bins)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        hist.T,
        origin="lower",
        aspect="auto",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="magma",
    )
    ax.set_xlabel("v_x")
    ax.set_ylabel("v_y")
    ax.set_title(f"Velocity mesh: {species} at ts={ts}")
    fig.colorbar(im, ax=ax, label="count")

    outname = f"velmesh_{species}_ts{ts}.png"
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
