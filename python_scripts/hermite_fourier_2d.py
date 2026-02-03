#!/usr/bin/env python3
"""
Hermite analysis in 2D: build f(x,vx) and f(y,vy), then plot Mx and My contours.

Usage:
  python hermite_fourier_2d.py <data_path> <species> [drift_vx] [drift_vy]
"""

import os
import sys
from os.path import join as pjoin
from multiprocessing import Pool, cpu_count

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp


DRIFT_VX_DEFAULT = 0.0
DRIFT_VY_DEFAULT = 0.0


def compute_hermite_basis(v, M):
    Nv = len(v)
    psi = np.zeros((M, Nv))
    psi[0] = np.pi ** (-0.25) * np.exp(-v ** 2 / 2)
    if M > 1:
        psi[1] = np.sqrt(2) * v * psi[0]
    for m in range(1, M - 1):
        psi[m + 1] = (
            np.sqrt(2 / (m + 1)) * v * psi[m]
            - np.sqrt(m / (m + 1)) * psi[m - 1]
        )
    return psi


def get_velocity_range(ts, path, species, drift_vx, drift_vy):
    with h5py.File(path, "r") as f:
        v = f[f"particle_{species}/vel{ts}"][:]
    vx = v[:, 0] - drift_vx
    vy = v[:, 1] - drift_vy
    return vx.min(), vx.max(), vy.min(), vy.max()


def process_timestep(args):
    idx, ts, path, species, x_edges, y_edges, v_edges_x, v_edges_y, psi_wx, psi_wy, dt_coeff, mfactor, drift_vx, drift_vy = args

    with h5py.File(path, "r") as f:
        pos = f[f"particle_{species}/pos{ts}"][:]
        vel = f[f"particle_{species}/vel{ts}"][:]

    x = pos[:, 0]
    y = pos[:, 1]
    vx = vel[:, 0] - drift_vx
    vy = vel[:, 1] - drift_vy

    f_xvx, _, _ = np.histogram2d(x, vx, bins=[x_edges, v_edges_x], density=True)
    f_yvy, _, _ = np.histogram2d(y, vy, bins=[y_edges, v_edges_y], density=True)

    f_mx = f_xvx @ psi_wx.T
    f_my = f_yvy @ psi_wy.T

    E_mx = np.mean(np.abs(f_mx) ** 2, axis=0)
    E_my = np.mean(np.abs(f_my) ** 2, axis=0)

    time_val = ts * dt_coeff * mfactor

    return idx, E_mx, E_my, time_val


def main():
    if len(sys.argv) < 3:
        print("Usage: python hermite_fourier_2d.py <data_path> <species> [drift_vx] [drift_vy]")
        sys.exit(1)

    data_path = sys.argv[1]
    species = sys.argv[2]
    drift_vx = float(sys.argv[3]) if len(sys.argv) > 3 else DRIFT_VX_DEFAULT
    drift_vy = float(sys.argv[4]) if len(sys.argv) > 4 else DRIFT_VY_DEFAULT

    path = pjoin(data_path, "result.h5")

    with h5py.File(path, "r") as f:
        meta = f["/metadata"].attrs
        nx = int(meta.get("Nx", 0))
        ny = int(meta.get("Ny", 0))
        NUM_TS = int(meta.get("NUM_TS", 0))
        write_int_phase = int(meta.get("write_int_phase", 1))
        DT_coeff = float(meta.get("DT_coeff", 1.0))
        wpe = float(meta.get("wpe", 1.0))
        wpi = float(meta.get("wpi", 1.0))
        normscheme = int(meta.get("norm_scheme", 1))

        group = f[f"particle_{species}"]
        vel_keys = sorted(int(k.replace("vel", "")) for k in group.keys() if k.startswith("vel"))

        if not vel_keys:
            raise RuntimeError(f"No velocity datasets found for species {species}")

        ts_list = vel_keys

        first_ts = ts_list[0]
        pos0 = group[f"pos{first_ts}"][:]
        x_min, x_max = pos0[:, 0].min(), pos0[:, 0].max()
        y_min, y_max = pos0[:, 1].min(), pos0[:, 1].max()

    mfactor = wpi / wpe if wpe != 0 else 1.0
    if normscheme in (1, 2):
        mfactor = 1.0

    if nx <= 1:
        nx = 64
    if ny <= 1:
        ny = 64

    x_edges = np.linspace(x_min, x_max, nx + 1)
    y_edges = np.linspace(y_min, y_max, ny + 1)

    # Scan velocity range
    ncores = max(1, min(cpu_count(), 6))
    with Pool(ncores) as p:
        vr = p.starmap(
            get_velocity_range,
            [(ts, path, species, drift_vx, drift_vy) for ts in ts_list],
        )

    vmin_x = min(v[0] for v in vr)
    vmax_x = max(v[1] for v in vr)
    vmin_y = min(v[2] for v in vr)
    vmax_y = max(v[3] for v in vr)

    v_range_x = vmax_x - vmin_x
    v_range_y = vmax_y - vmin_y

    vmin_x -= 0.1 * v_range_x
    vmax_x += 0.1 * v_range_x
    vmin_y -= 0.1 * v_range_y
    vmax_y += 0.1 * v_range_y

    Nv = 1024
    v_edges_x = np.linspace(vmin_x, vmax_x, Nv)
    v_edges_y = np.linspace(vmin_y, vmax_y, Nv)

    v_centers_x = 0.5 * (v_edges_x[:-1] + v_edges_x[1:])
    v_centers_y = 0.5 * (v_edges_y[:-1] + v_edges_y[1:])

    dvx = v_centers_x[1] - v_centers_x[0]
    dvy = v_centers_y[1] - v_centers_y[0]

    M_max = 250

    psi_x = compute_hermite_basis(v_centers_x, M_max)
    psi_y = compute_hermite_basis(v_centers_y, M_max)

    psi_wx = psi_x * dvx
    psi_wy = psi_y * dvy

    DATA_TS = len(ts_list)
    E_mx_time = np.zeros((DATA_TS, M_max))
    E_my_time = np.zeros((DATA_TS, M_max))
    time_array = np.zeros(DATA_TS)

    work = [
        (
            i,
            ts,
            path,
            species,
            x_edges,
            y_edges,
            v_edges_x,
            v_edges_y,
            psi_wx,
            psi_wy,
            DT_coeff,
            mfactor,
            drift_vx,
            drift_vy,
        )
        for i, ts in enumerate(ts_list)
    ]

    with Pool(ncores) as p:
        for idx, E_mx, E_my, tt in p.imap_unordered(process_timestep, work):
            E_mx_time[idx] = E_mx
            E_my_time[idx] = E_my
            time_array[idx] = tt

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
    plt.contourf(time_array, np.arange(M_max), np.log10(E_mx_time.T + 1e-20), 40, cmap="RdYlBu_r")
    plt.xlabel(r"$\omega_{pi} t$")
    plt.ylabel(r"$m_x$")
    plt.colorbar(label=r"$\log|f_{m_x}|^2$")
    plt.savefig(pjoin(outdir, "hermite_time_mx.png"), dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
    plt.contourf(time_array, np.arange(M_max), np.log10(E_my_time.T + 1e-20), 40, cmap="RdYlBu_r")
    plt.xlabel(r"$\omega_{pi} t$")
    plt.ylabel(r"$m_y$")
    plt.colorbar(label=r"$\log|f_{m_y}|^2$")
    plt.savefig(pjoin(outdir, "hermite_time_my.png"), dpi=dpi)
    plt.close()

    print(f"Saved plots to {outdir}")


if __name__ == "__main__":
    main()
