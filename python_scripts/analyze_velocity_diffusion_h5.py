#!/usr/bin/env python3
"""
Plot velocity diffusion diagnostics from a saved diffusion HDF5 file.

Usage:
  python analyze_velocity_diffusion_h5.py <path_to_simulation> <particle_type> [component] [V0]

Arguments:
  component : x | y | mag | both (default: x)
  V0        : velocity bin center to focus on (default: v_phase or mid-bin)
"""

import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp


def parse_component(comp_raw):
    comp = comp_raw.lower()
    if comp in ("x", "vx"):
        return "x"
    if comp in ("y", "vy"):
        return "y"
    if comp in ("mag", "m", "speed"):
        return "mag"
    if comp in ("both", "xy", "x+y"):
        return "both"
    raise ValueError("component must be x, y, mag, or both")


def trapz(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def load_diffusion(fname):
    with h5py.File(fname, "r") as f:
        axes = f["axes"]
        meta = f["metadata"].attrs

        v_centers = axes["v_centers"][:]
        tau_phys = axes["tau_phys"][:]
        tau_steps = axes["tau_steps"][:]

        all_pdfs = f["pdfs/P_dv_tau_v"][:]
        sigma = f["gaussian_fit/sigma"][:]
        mu = f["gaussian_fit/mu"][:]

        D_var = f["local_diffusion/D_from_var"][:]
        D_msd = f["local_diffusion/D_from_msd"][:]
        D_sigma = f["local_diffusion/D_from_sigma"][:]

        pdf_centers_grp = axes["pdf_centers"]
        pdf_centers = {k: pdf_centers_grp[k][:] for k in pdf_centers_grp.keys()}

        v_phase = float(meta.get("v_phase", 0.0))
        Nfit = int(meta.get("Nfit", min(6, len(tau_phys))))

        all_var = f["statistics/variance"][:]
        all_sigma_fit = f["gaussian_fit/sigma"][:]

    return {
        "v_centers": v_centers,
        "tau_phys": tau_phys,
        "tau_steps": tau_steps,
        "all_pdfs": all_pdfs,
        "sigma": sigma,
        "mu": mu,
        "D_var": D_var,
        "D_msd": D_msd,
        "D_sigma": D_sigma,
        "pdf_centers": pdf_centers,
        "v_phase": v_phase,
        "Nfit": Nfit,
        "all_var": all_var,
        "all_sigma_fit": all_sigma_fit,
    }


def compute_kurtosis(all_pdfs, pdf_centers):
    nv_bins, n_tau, _ = all_pdfs.shape
    kurt = np.full((nv_bins, n_tau), np.nan)

    for ib in range(nv_bins):
        dv = pdf_centers.get(f"bin_{ib}")
        if dv is None or np.all(~np.isfinite(dv)):
            continue
        for it in range(n_tau):
            P = all_pdfs[ib, it]
            norm = trapz(P, dv)
            if norm <= 0:
                continue
            P = P / norm
            mean = trapz(dv * P, dv)
            var = trapz((dv - mean) ** 2 * P, dv)
            if var <= 0:
                continue
            kurt[ib, it] = trapz((dv - mean) ** 4 * P, dv) / var ** 2

    return kurt


if len(sys.argv) < 3:
    print("Usage: python analyze_velocity_diffusion_h5.py <path_to_simulation> <particle_type> [component] [V0]")
    sys.exit(1)

path = sys.argv[1]
ptype = sys.argv[2]
component = parse_component(sys.argv[3]) if len(sys.argv) > 3 else "x"

analysis_dir = pjoin(path, "new_analysis")

path_fig = pjoin(path, "diffusion_analysis_plots")
os.makedirs(path_fig, exist_ok=True)

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

if component == "both":
    data = {}
    for comp in ("x", "y"):
        fname = pjoin(analysis_dir, f"local_velocity_diffusion_{ptype}_{comp}.h5")
        data[comp] = load_diffusion(fname)

    v_centers = data["x"]["v_centers"]
    tau_phys = data["x"]["tau_phys"]

    sigma_map_x = data["x"]["all_var"]
    sigma_map_y = data["y"]["all_var"]

    kurt_x = compute_kurtosis(data["x"]["all_pdfs"], data["x"]["pdf_centers"])
    kurt_y = compute_kurtosis(data["y"]["all_pdfs"], data["y"]["pdf_centers"])

    fig, axes = plt.subplots(1, 2, figsize=figsize / 20.0, constrained_layout=True, dpi=ppi)
    im0 = axes[0].imshow(
        sigma_map_x,
        aspect="auto",
        origin="lower",
        interpolation="bilinear",
        extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
        cmap="RdBu_r",
    )
    axes[0].set_title("sigma^2 map (vx)")
    axes[0].set_xlabel(r"$\tau$")
    axes[0].set_ylabel(r"$v$")

    im1 = axes[1].imshow(
        sigma_map_y,
        aspect="auto",
        origin="lower",
        interpolation="bilinear",
        extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
        cmap="RdBu_r",
    )
    axes[1].set_title("sigma^2 map (vy)")
    axes[1].set_xlabel(r"$\tau$")
    axes[1].set_ylabel(r"$v$")

    fig.colorbar(im1, ax=axes.ravel().tolist(), label=r"$\sigma^2$")
    plt.savefig(pjoin(path_fig, "sigma_map_xy.pdf"), dpi=dpi)
    plt.show()

    fig, axes = plt.subplots(1, 2, figsize=figsize / 20.0, constrained_layout=True, dpi=ppi)
    im0 = axes[0].imshow(
        kurt_x,
        aspect="auto",
        origin="lower",
        interpolation="bilinear",
        extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
    )
    axes[0].set_title("kurtosis map (vx)")
    axes[0].set_xlabel(r"$\tau$")
    axes[0].set_ylabel(r"$v$")

    im1 = axes[1].imshow(
        kurt_y,
        aspect="auto",
        origin="lower",
        interpolation="bilinear",
        extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
    )
    axes[1].set_title("kurtosis map (vy)")
    axes[1].set_xlabel(r"$\tau$")
    axes[1].set_ylabel(r"$v$")

    fig.colorbar(im1, ax=axes.ravel().tolist(), label="Kurtosis k")
    plt.savefig(pjoin(path_fig, "kurtosis_map_xy.pdf"), dpi=dpi)
    plt.show()

    sys.exit(0)

fname = pjoin(analysis_dir, f"local_velocity_diffusion_{ptype}_{component}.h5")
data = load_diffusion(fname)

v_centers = data["v_centers"]
tau_phys = data["tau_phys"]
tau_steps = data["tau_steps"]

all_pdfs = data["all_pdfs"]
sigma = data["sigma"]
mu = data["mu"]

D_var = data["D_var"]
D_msd = data["D_msd"]
D_sigma = data["D_sigma"]

pdf_centers = data["pdf_centers"]
v_phase = data["v_phase"]
Nfit = data["Nfit"]
all_var = data["all_var"]
all_sigma_fit = data["all_sigma_fit"]

nv_bins, n_tau, _ = all_pdfs.shape

if len(sys.argv) > 4:
    V0 = float(sys.argv[4])
elif v_phase != 0.0:
    V0 = v_phase
else:
    V0 = v_centers[len(v_centers) // 2]

# Diffusion coefficients vs velocity
plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
plt.plot(v_centers, D_msd, "-o", label=r"$D_{1}(v)$")
plt.plot(v_centers, D_sigma, "-s", label=r"$D_{2}(v)$")
plt.xlim([v_centers[0], v_centers[-1]])
plt.xlabel(r"$v$")
plt.ylabel(r"$D(v)$")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_fig, f"diffusion_coefficient_vs_velocity_{component}.pdf"), dpi=dpi)
plt.show()

# Variance vs tau for all bins
plt.figure()
for ib in range(0, nv_bins, 4):
    plt.plot(tau_phys, all_var[ib], alpha=0.4)
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\sigma^2(\Delta v)$")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# Diffusion exponent n(v): sigma^2 ~ tau^n
n_v = np.zeros(nv_bins)

for ib in range(nv_bins):
    y = sigma[ib, :Nfit] ** 2
    mask = np.isfinite(y)
    if np.count_nonzero(mask) >= 2:
        n_v[ib] = np.polyfit(np.log(tau_phys[:Nfit][mask]), np.log(y[mask]), 1)[0]
    else:
        n_v[ib] = np.nan

plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
plt.plot(v_centers, n_v, "--", color="blue")
plt.scatter(v_centers, n_v, color="blue", s=40)
plt.xlabel(r"$v$")
plt.ylabel(r"$n$  ($\sigma^2 \sim \tau^{\,n}$)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(pjoin(path_fig, f"diffusion_exponent_vs_velocity_{component}.pdf"), dpi=dpi)
plt.show()

# Non-Gaussianity: kurtosis
kurt = compute_kurtosis(all_pdfs, pdf_centers)

plt.figure()
plt.imshow(
    kurt,
    aspect="auto",
    origin="lower",
    interpolation="bilinear",
    extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
)
plt.colorbar(label="Kurtosis k")
plt.xlabel(r"$\tau$")
plt.ylabel(r"$v$")
plt.title("Non-Gaussianity map")
plt.tight_layout()
plt.show()

# Select velocity bin closest to V0
ib = int(np.argmin(np.abs(v_centers - V0)))
print("bin vel close to V0", v_centers[ib])

dv = pdf_centers.get(f"bin_{ib}")
if dv is not None and np.any(np.isfinite(dv)):
    tau_list = [1, 4, 8, 12, 16, 20]
    tau_list = [t for t in tau_list if t < n_tau]

    plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)

    for it in tau_list:
        P = all_pdfs[ib, it]
        norm = trapz(P, dv)
        if norm > 0:
            P = P / norm
        plt.plot(dv, P, lw=1.8, label=rf"$\tau = {it:.2f}$")

    plt.xlabel(r"$\Delta v$")
    plt.ylabel(r"$P(\Delta v,\tau)$")
    plt.yscale("log")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(pjoin(path_fig, f"PDF_vs_tau_{component}.pdf"), dpi=dpi)
    plt.show()

plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
plt.imshow(
    all_var,
    aspect="auto",
    origin="lower",
    interpolation="bilinear",
    extent=[tau_phys[0], tau_phys[-1], v_centers[0], v_centers[-1]],
    cmap="RdBu_r",
)
plt.colorbar(label=r"$\sigma^2$")
plt.xlabel(r"$\tau$")
plt.ylabel(r"$v$")
plt.savefig(pjoin(path_fig, f"sigma_map_{component}.pdf"), dpi=dpi)
plt.show()

# Sigma^2 vs tau for a single bin
msd_array = all_var[ib]
sigma2 = all_sigma_fit[ib] ** 2

mask_msd = np.isfinite(msd_array)
mask_sig = np.isfinite(sigma2)

if np.count_nonzero(mask_msd) >= 2 and np.count_nonzero(mask_sig) >= 2:
    slope_msd, intercept_msd = np.polyfit(tau_phys[:Nfit], msd_array[:Nfit], 1)
    slope_sigma, intercept_sigma = np.polyfit(tau_phys[:Nfit], sigma2[:Nfit], 1)

    sigma2_msd_fit = slope_msd * tau_phys[:Nfit] + intercept_msd
    sigma2_gauss_fit = slope_sigma * tau_phys[:Nfit] + intercept_sigma

    plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)
    plt.scatter(tau_phys[:Nfit], msd_array[:Nfit], s=30, color="blue", label=r"$\sigma^2_{1}$")
    plt.scatter(tau_phys[:Nfit], sigma2[:Nfit], s=30, color="red", label=r"$\sigma^2_{2}$")
    plt.plot(tau_phys[:Nfit], sigma2_msd_fit, color="blue", linewidth=2)
    plt.plot(tau_phys[:Nfit], sigma2_gauss_fit, color="red", linewidth=2)
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$\sigma^2$")
    plt.grid(alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(pjoin(path_fig, f"sigma2_vs_tau_bin{ib}_{component}.pdf"), dpi=900)
    plt.show()

# PDF self-similarity test
if dv is not None and np.any(np.isfinite(dv)):
    tau_list = [1, 4, 8, 12, 16, 20]
    tau_list = [t for t in tau_list if t < n_tau]

    plt.figure(figsize=figsize / 25.4, constrained_layout=True, dpi=ppi)

    for it in tau_list:
        P = all_pdfs[ib, it]
        s = sigma[ib, it]
        if not np.isfinite(s) or s == 0:
            continue
        plt.plot(dv / s, P * s, lw=1.8, alpha=0.7, label=rf"$\tau = {it:.2f}$")

    x = np.linspace(-5, 5, 400)
    plt.plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), "k--", lw=2, label="Gaussian")
    plt.xlabel(r"$\Delta v / \sigma$")
    plt.ylabel(r"$\sigma\, P(\Delta v,\tau)$")
    plt.yscale("log")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(pjoin(path_fig, f"PDF_self_similarity_test_{component}.pdf"), dpi=dpi)
    plt.show()
