import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 particle_trajectory.py <path> <species> <particle_index> [mode]")
    print("mode: xy | x | y (default: xy)")


if len(sys.argv) < 4:
    usage()
    sys.exit(1)

path = sys.argv[1]
species = sys.argv[2]
particle_index = int(sys.argv[3])
mode = sys.argv[4].lower() if len(sys.argv) >= 5 else "xy"

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

    timesteps = list(range(0, num_ts + 1, write_int_phase))
    xs, ys, ts = [], [], []

    for t in timesteps:
        ds_name = f"{phase_group}/pos{t}"
        if ds_name not in f:
            continue
        pos = f[ds_name][:]
        if particle_index >= pos.shape[0]:
            print(f"Particle index {particle_index} out of range at ts={t}")
            sys.exit(1)
        xs.append(pos[particle_index, 0])
        ys.append(pos[particle_index, 1])
        ts.append(t)

    if not ts:
        print("No particle data found")
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(8, 6))
    if mode == "xy":
        ax.plot(xs, ys, marker="o", markersize=3, linewidth=1.2)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"Trajectory (x-y) of particle {particle_index} ({species})")
        outname = f"trajectory_xy_{species}_{particle_index}.png"
    elif mode == "x":
        ax.plot(ts, xs, marker="o", markersize=3, linewidth=1.2)
        ax.set_xlabel("Time step")
        ax.set_ylabel("x")
        ax.set_title(f"x(t) of particle {particle_index} ({species})")
        outname = f"trajectory_x_{species}_{particle_index}.png"
    elif mode == "y":
        ax.plot(ts, ys, marker="o", markersize=3, linewidth=1.2)
        ax.set_xlabel("Time step")
        ax.set_ylabel("y")
        ax.set_title(f"y(t) of particle {particle_index} ({species})")
        outname = f"trajectory_y_{species}_{particle_index}.png"
    else:
        print(f"Unknown mode: {mode}")
        usage()
        sys.exit(1)

    ax.grid(True, alpha=0.3)
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
