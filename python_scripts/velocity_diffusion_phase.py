#!/usr/bin/env python3
"""
Compute local velocity diffusion from phase-space data only (no field analysis).

Usage:
  python velocity_diffusion_phase.py <path_to_simulation> <particle_type> [component] [t_start] [t_end] [nv_bins]

Arguments:
  component : x | y | mag (default: x)
  t_start   : snapshot index start (default: 0)
  t_end     : snapshot index end (exclusive, default: last)
  nv_bins   : number of velocity bins (default: 50)
"""

import os
import sys
from os.path import join as pjoin

import h5py
import numpy as np
from scipy.optimize import curve_fit


def gauss(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def parse_component(comp_raw):
    comp = comp_raw.lower()
    if comp in ("x", "vx"):
        return "x"
    if comp in ("y", "vy"):
        return "y"
    if comp in ("mag", "m", "speed"):
        return "mag"
    raise ValueError("component must be x, y, or mag")


def extract_vel_component(v, comp):
    if comp == "x":
        return v[:, 0]
    if comp == "y":
        return v[:, 1]
    return np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)


if len(sys.argv) < 3:
    print("Usage: python velocity_diffusion_phase.py <path_to_simulation> <particle_type> [component] [t_start] [t_end] [nv_bins]")
    sys.exit(1)

path = sys.argv[1]
particle_type = sys.argv[2]
component = parse_component(sys.argv[3]) if len(sys.argv) > 3 else "x"

t_start = int(sys.argv[4]) if len(sys.argv) > 4 else 0
t_end = int(sys.argv[5]) if len(sys.argv) > 5 else -1
nv_bins = int(sys.argv[6]) if len(sys.argv) > 6 else 50

file_name = "result.h5"
analysis_dir = pjoin(path, "new_analysis")
os.makedirs(analysis_dir, exist_ok=True)

with h5py.File(pjoin(path, file_name), "r") as f:
    meta = f["/metadata"].attrs

    NUM_TS = int(meta.get("NUM_TS", 0))
    write_interval_phase = int(meta.get("write_int_phase", 1))
    DT_coeff = float(meta.get("DT_coeff", 1.0))
    wpe = float(meta.get("wpe", 1.0))
    wpi = float(meta.get("wpi", 1.0))
    normscheme = int(meta.get("norm_scheme", 1))

    mfactor = wpi / wpe if wpe != 0 else 1.0
    if normscheme in (1, 2):
        mfactor = 1.0

    dt_phase = write_interval_phase * DT_coeff * mfactor

    particle_group_name = f"particle_{particle_type}"
    if particle_group_name not in f:
        raise RuntimeError(f"Group {particle_group_name} not found in HDF5 file")

    gpart = f[particle_group_name]
    vel_keys = [k for k in gpart.keys() if k.startswith("vel")]
    if not vel_keys:
        raise RuntimeError(f"No velocity datasets found in {particle_group_name}")

    ts_list = sorted(int(k.replace("vel", "")) for k in vel_keys)

    vel_list = []
    n_list = []
    for ts in ts_list:
        v = gpart[f"vel{ts}"][:]
        v_comp = extract_vel_component(v, component).astype(float)
        vel_list.append(v_comp)
        n_list.append(v_comp.size)

    min_n = min(n_list)
    max_n = max(n_list)
    if min_n != max_n:
        print(f"Warning: particle count changes across snapshots (min={min_n}, max={max_n}). Truncating to min.")

    vel = np.stack([v[:min_n] for v in vel_list], axis=0)
    nt, npart = vel.shape

    if t_end < 0 or t_end > nt:
        t_end = nt
    t_start = max(t_start, 0)
    if t_end <= t_start:
        raise RuntimeError("t_end must be greater than t_start")

    ti_idxs = np.arange(t_start, t_end)

    v_flat = vel[ti_idxs].flatten()
    vmin = np.percentile(v_flat, 5)
    vmax = np.percentile(v_flat, 95)
    if vmin == vmax:
        raise RuntimeError("Velocity range is zero; cannot bin.")

    v_edges = np.linspace(vmin, vmax, nv_bins + 1)
    v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])

    part1 = np.arange(1, 6)
    part2 = np.arange(6, 31, 2)
    part3 = np.arange(32, 101, 8)
    tau_steps = np.unique(np.concatenate([part1, part2, part3]))

    max_tau = t_end - t_start - 1
    tau_steps = tau_steps[tau_steps <= max_tau]
    if tau_steps.size == 0:
        raise RuntimeError("Not enough snapshots to compute any tau steps.")

    tau_phys = tau_steps * dt_phase

    n_pdf_bins = 100

    all_pdfs = np.zeros((nv_bins, len(tau_steps), n_pdf_bins))
    all_sigma_fit = np.full((nv_bins, len(tau_steps)), np.nan)
    all_mu_fit = np.full((nv_bins, len(tau_steps)), np.nan)
    all_amp_fit = np.full((nv_bins, len(tau_steps)), np.nan)
    all_D2_from_sigma_tau = np.full((nv_bins, len(tau_steps)), np.nan)

    all_var = np.full((nv_bins, len(tau_steps)), np.nan)
    all_sigma_msd = np.full((nv_bins, len(tau_steps)), np.nan)

    D1 = np.full(nv_bins, np.nan)
    D2 = np.full(nv_bins, np.nan)
    D2_from_sigma = np.full(nv_bins, np.nan)

    pdf_centers_storage = []

    Nfit = max(len(tau_steps) // 4, 6)
    Nfit = min(Nfit, len(tau_steps))

    for ib in range(nv_bins):
        v_low = v_edges[ib]
        v_high = v_edges[ib + 1]

        dv_per_tau = {tau: [] for tau in tau_steps}

        for ti in ti_idxs:
            v_t = vel[ti]
            in_bin_mask = (v_t >= v_low) & (v_t < v_high)
            idxs = np.flatnonzero(in_bin_mask)
            if idxs.size == 0:
                continue

            v_initial = v_t[idxs]
            for tau in tau_steps:
                t2 = ti + tau
                if t2 >= nt:
                    continue
                v_final = vel[t2, idxs]
                dv_per_tau[tau].append(v_final - v_initial)

        all_dv_samples_for_bin = []
        for tau in tau_steps:
            if dv_per_tau[tau]:
                all_dv_samples_for_bin.append(np.concatenate(dv_per_tau[tau]))

        if not all_dv_samples_for_bin:
            pdf_centers_storage.append(np.full(n_pdf_bins, np.nan))
            continue

        all_dv_concat = np.concatenate(all_dv_samples_for_bin)
        dv_min = np.min(all_dv_concat)
        dv_max = np.max(all_dv_concat)
        dv_range = dv_max - dv_min
        if dv_range == 0:
            pdf_centers_storage.append(np.full(n_pdf_bins, np.nan))
            continue

        buffer_width = 0.05 * dv_range
        pdf_bins = np.linspace(dv_min - buffer_width, dv_max + buffer_width, n_pdf_bins + 1)
        pdf_bin_centers = 0.5 * (pdf_bins[:-1] + pdf_bins[1:])
        pdf_centers_storage.append(pdf_bin_centers)

        for i_tau, tau in enumerate(tau_steps):
            dv_list = dv_per_tau[tau]
            if not dv_list:
                continue

            dv_cat = np.concatenate(dv_list)
            msd_val = np.mean(dv_cat ** 2)
            var_val = np.var(dv_cat)

            all_var[ib, i_tau] = var_val
            all_sigma_msd[ib, i_tau] = np.sqrt(msd_val)

            hist_vals, _ = np.histogram(dv_cat, bins=pdf_bins, density=True)
            all_pdfs[ib, i_tau, :] = hist_vals

            peak_idx = int(np.argmax(hist_vals))
            mu0 = pdf_bin_centers[peak_idx]
            a0 = hist_vals[peak_idx]
            sigma0 = np.std(dv_cat) if np.std(dv_cat) > 0 else 1e-12

            p0 = [a0, mu0, sigma0]
            bounds = ([0.0, pdf_bins[0], 1e-12], [np.inf, pdf_bins[-1], np.inf])

            try:
                popt, _ = curve_fit(gauss, pdf_bin_centers, hist_vals, p0=p0, bounds=bounds, maxfev=8000)
                a_fit, mu_f, sigma_f = popt
                all_amp_fit[ib, i_tau] = a_fit
                all_mu_fit[ib, i_tau] = mu_f
                all_sigma_fit[ib, i_tau] = abs(sigma_f)
                tau_val = tau_phys[i_tau]
                if tau_val != 0:
                    all_D2_from_sigma_tau[ib, i_tau] = (sigma_f ** 2) / (2.0 * tau_val)
            except RuntimeError:
                continue

        # Diffusion coefficients for this bin
        tau_fit = tau_phys[:Nfit]

        var_slice = all_var[ib, :Nfit]
        msd_slice = all_sigma_msd[ib, :Nfit] ** 2
        sigma_fit_slice = all_sigma_fit[ib, :Nfit]

        mask_var = np.isfinite(var_slice)
        if np.count_nonzero(mask_var) >= 2:
            slope_var = np.polyfit(tau_fit[mask_var], var_slice[mask_var], 1)[0]
            D1[ib] = slope_var / 2.0

        mask_msd = np.isfinite(msd_slice)
        if np.count_nonzero(mask_msd) >= 2:
            slope_msd = np.polyfit(tau_fit[mask_msd], msd_slice[mask_msd], 1)[0]
            D2[ib] = slope_msd / 2.0

        sigma2 = sigma_fit_slice ** 2
        mask_sigma = np.isfinite(sigma2)
        if np.count_nonzero(mask_sigma) >= 2:
            slope_sigma = np.polyfit(tau_fit[mask_sigma], sigma2[mask_sigma], 1)[0]
            D2_from_sigma[ib] = slope_sigma / 2.0

    # Save HDF5 output
    save_h5 = pjoin(analysis_dir, f"local_velocity_diffusion_{particle_type}_{component}.h5")
    with h5py.File(save_h5, "w") as hf:
        gmeta = hf.create_group("metadata")
        gmeta.attrs["particle_type"] = particle_type
        gmeta.attrs["component"] = component
        gmeta.attrs["nv_bins"] = nv_bins
        gmeta.attrs["n_tau"] = len(tau_steps)
        gmeta.attrs["dt_phase"] = dt_phase
        gmeta.attrs["v_phase"] = 0.0
        gmeta.attrs["analysis_window_start"] = t_start
        gmeta.attrs["analysis_window_end"] = t_end
        gmeta.attrs["Nfit"] = Nfit
        gmeta.attrs["npart_min"] = min_n
        gmeta.attrs["npart_max"] = max_n

        axes = hf.create_group("axes")
        axes.create_dataset("v_edges", data=v_edges)
        axes.create_dataset("v_centers", data=v_centers)
        axes.create_dataset("tau_steps", data=tau_steps)
        axes.create_dataset("tau_phys", data=tau_phys)

        pdf_axes = axes.create_group("pdf_centers")
        for ib in range(nv_bins):
            pdf_axes.create_dataset(f"bin_{ib}", data=pdf_centers_storage[ib])

        gpdf = hf.create_group("pdfs")
        gpdf.create_dataset("P_dv_tau_v", data=all_pdfs, compression="gzip")

        gfit = hf.create_group("gaussian_fit")
        gfit.create_dataset("A", data=all_amp_fit)
        gfit.create_dataset("mu", data=all_mu_fit)
        gfit.create_dataset("sigma", data=all_sigma_fit)

        gDsig = hf.create_group("diffusion_from_sigma_tau")
        gDsig.create_dataset("D_tau", data=all_D2_from_sigma_tau)

        gD = hf.create_group("local_diffusion")
        gD.create_dataset("D_from_var", data=D1)
        gD.create_dataset("D_from_msd", data=D2)
        gD.create_dataset("D_from_sigma", data=D2_from_sigma)

        gstats = hf.create_group("statistics")
        gstats.create_dataset("variance", data=all_var)
        gstats.create_dataset("sigma_from_msd", data=all_sigma_msd)

    print(f"Saved diffusion analysis to: {save_h5}")
