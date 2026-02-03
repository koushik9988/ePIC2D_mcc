import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt


def usage():
    print("Usage: python3 potential_energy.py <path>")


if len(sys.argv) < 2:
    usage()
    sys.exit(1)

path = sys.argv[1]
file_name = "result.h5"
plot_dir = pjoin(path, "plots")
os.makedirs(plot_dir, exist_ok=True)

with h5py.File(pjoin(path, file_name), "r") as f:
    spno = f["/metadata"].attrs.get("spno", None)
    if spno is None:
        print("Missing /metadata spno")
        sys.exit(1)

    data = f["time_var/kinetic_energy"]
    ts = data[:, 0]

    pex = data[:, 3 * spno + 1]
    pey = data[:, 3 * spno + 2]
    pe = pex + pey

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ts, pex, label="PEx")
    ax.plot(ts, pey, label="PEy")
    ax.plot(ts, pe, label="PE total", color="black")
    ax.set_xlabel("Time")
    ax.set_ylabel("Potential Energy")
    ax.set_title("Potential Energy vs Time")
    ax.legend()
    ax.grid(True, alpha=0.3)

    outname = "potential_energy.png"
    fig.savefig(pjoin(plot_dir, outname), dpi=300)
    print(f"Saved plot to {pjoin(plot_dir, outname)}")
    plt.show()
