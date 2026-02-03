#!/usr/bin/env python3
"""
2D FFT of f(vx, vy, t) using sequential 1D FFTs (vx then vy) to limit memory.
Averages spectrum over time and plots Mx-My contour.

Usage:
  python fft_vx_vy_spectrum.py <data_path> <species> [nv_bins] [t_start] [t_end] [stride] [drift_vx] [drift_vy]

Defaults:
  nv_bins=256, t_start=0, t_end=-1 (all), stride=1, drift_vx=0, drift_vy=0
"""

import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp


def parse_int(arg, default):
    if arg is None:
        return default
    return int(arg)


def parse_float(arg, default):
    if arg is None:
        return default
    return float(arg)


def main():
    if len(sys.argv) < 3:
        print("Usage: python fft_vx_vy_spectrum.py <data_path> <species> [nv_bins] [t_start] [t_end] [stride] [drift_vx] [drift_vy]")
        sys.exit(1)

    data_path = sys.argv[1]
    species = sys.argv[2]

    nv_bins = parse_int(sys.argv[3] if len(sys.argv) > 3 else None, 256)
    t_start = parse_int(sys.argv[4] if len(sys.argv) > 4 else None, 0)
    t_end = parse_int(sys.argv[5] if len(sys.argv) > 5 else None, -1)
    stride = parse_int(sys.argv[6] if len(sys.argv) > 6 else None, 1)
    drift_vx = parse_float(sys.argv[7] if len(sys.argv) > 7 else None, 0.0)
    drift_vy = parse_float(sys.argv[8] if len(sys.argv) > 8 else None, 0.0)

    path = pjoin(data_path, "result.h5")

    with h5py.File(path, "r") as f:
        group_name = f"particle_{species}"
        if group_name not in f:
            raise RuntimeError(f"Group {group_name} not found")
        g = f[group_name]
        vel_keys = sorted(int(k.replace("vel", "")) for k in g.keys() if k.startswith("vel"))
        if not vel_keys:
            raise RuntimeError(f"No velocity datasets for species {species}")

        if t_end < 0 or t_end > len(vel_keys):
            t_end = len(vel_keys)
        t_start = max(0, t_start)

        ts_list = vel_keys[t_start:t_end:stride]

        # Scan velocity ranges
        vmin_x = np.inf
        vmax_x = -np.inf
        vmin_y = np.inf
        vmax_y = -np.inf

        for ts in ts_list:
            vel = g[f"vel{ts}"][:]
            vx = vel[:, 0] - drift_vx
            vy = vel[:, 1] - drift_vy
            vmin_x = min(vmin_x, vx.min())
            vmax_x = max(vmax_x, vx.max())
            vmin_y = min(vmin_y, vy.min())
            vmax_y = max(vmax_y, vy.max())

    # Add buffer
    vx_range = vmax_x - vmin_x
    vy_range = vmax_y - vmin_y
    vmin_x -= 0.1 * vx_range
    vmax_x += 0.1 * vx_range
    vmin_y -= 0.1 * vy_range
    vmax_y += 0.1 * vy_range

    v_edges_x = np.linspace(vmin_x, vmax_x, nv_bins + 1)
    v_edges_y = np.linspace(vmin_y, vmax_y, nv_bins + 1)

    dvx = v_edges_x[1] - v_edges_x[0]
    dvy = v_edges_y[1] - v_edges_y[0]

    nkx = nv_bins // 2 + 1
    nky = nv_bins

    spec_sum = np.zeros((nkx, nky), dtype=np.float64)
    count = 0

    with h5py.File(path, "r") as f:
        g = f[f"particle_{species}"]
        for ts in ts_list:
            vel = g[f"vel{ts}"][:]
            vx = vel[:, 0] - drift_vx
            vy = vel[:, 1] - drift_vy

            f_vxvy, _, _ = np.histogram2d(vx, vy, bins=[v_edges_x, v_edges_y], density=True)

            # Sequential FFT: vx first, then vy
            F_vx = np.fft.rfft(f_vxvy, axis=0)
            F_vxvy = np.fft.fft(F_vx, axis=1)

            power = np.abs(F_vxvy) ** 2
            spec_sum += power
            count += 1

    if count == 0:
        raise RuntimeError("No snapshots processed. Check t_start/t_end/stride.")

    spec_avg = spec_sum / count

    mx = np.arange(spec_avg.shape[0])
    my = np.arange(spec_avg.shape[1])

    outdir = pjoin(data_path, "paper_plots")
    os.makedirs(outdir, exist_ok=True)

    figsize = np.array([150, 150 / 1.618])
    dpi = 300
    ppi = np.sqrt(1920**2 + 1200**2) / 24

    mp.rc("text", usetex=False)
    mp.rc("font", family="sans-serif", size=10, serif="Computer Modern Roman")
    mp.rc("axes", titlesize=10)
    mp.rc("axes", labelsize=10)
    mp.rc("xtick", labelsize=10)
    mp.rc("ytick", labelsize=10)
    mp.rc("legend", fontsize=10)

    plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
    plt.contourf(mx, my, np.log10(spec_avg.T + 1e-20), 40, cmap="RdYlBu_r")
    plt.xlabel(r"$m_x$")
    plt.ylabel(r"$m_y$")
    plt.colorbar(label=r"$\log |F(m_x,m_y)|^2$")
    plt.savefig(pjoin(outdir, f"fft_vx_vy_spectrum_{species}.png"), dpi=dpi)
    plt.close()

    save_h5 = pjoin(data_path, "new_analysis", f"fft_vx_vy_{species}.h5")
    os.makedirs(pjoin(data_path, "new_analysis"), exist_ok=True)
    with h5py.File(save_h5, "w") as hf:
        hf.create_dataset("power", data=spec_avg)
        hf.create_dataset("mx", data=mx)
        hf.create_dataset("my", data=my)
        hf.create_dataset("v_edges_x", data=v_edges_x)
        hf.create_dataset("v_edges_y", data=v_edges_y)
        hf.attrs["nv_bins"] = nv_bins
        hf.attrs["t_start"] = t_start
        hf.attrs["t_end"] = t_end
        hf.attrs["stride"] = stride
        hf.attrs["drift_vx"] = drift_vx
        hf.attrs["drift_vy"] = drift_vy
        hf.attrs["dvx"] = dvx
        hf.attrs["dvy"] = dvy
        hf.attrs["count"] = count

    print(f"Saved contour plot to {outdir}")
    print(f"Saved spectrum to {save_h5}")


if __name__ == "__main__":
    main()
