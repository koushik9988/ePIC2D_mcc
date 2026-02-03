import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from os.path import join as pjoin
from scipy.signal import find_peaks
from scipy.constants import value as constants

# === Input Arguments ===
if len(sys.argv) != 3:
    print("Usage: python script.py <path_to_data> <slice_index>")
    sys.exit(1)

path = sys.argv[1]
slice_index = int(sys.argv[2])

file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')

# ----------------------------- Read Metadata -----------------------------
metadata = f['/metadata']
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']

# Constants
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
AMU = constants('atomic mass constant')

# Metadata
nx = metadata.attrs['Nx']
ny = metadata.attrs['Ny']
NUM_TS = metadata.attrs['NUM_TS']
write_interval = metadata.attrs['write_int']
DT_coeff = metadata.attrs['DT_coeff']
wpe = metadata.attrs['wpe']
LDe = metadata.attrs['LDe']
wpi = metadata.attrs['wpi']
Te = metadata_electron.attrs['temperature'] * e
mi = metadata_ion.attrs['mass']
Bz = 0.0001  # Assumed magnetic field

DT = DT_coeff * (1.0 / wpe)

# Validate slice_index
if slice_index >= ny or slice_index >= nx:
    print(f"Error: slice_index {slice_index} out of bounds for nx={nx}, ny={ny}")
    sys.exit(1)

# ----------------------------- Load Electric Field Data -----------------------------
def load_field_data(field_component, slice_index):
    field_data = []
    field_group = f[f'fielddata/{field_component}']
    time_steps = sorted(map(int, field_group.keys()))
    for time_step in time_steps:
        EF_flat = field_group[str(time_step)][:]
        EF_2D = EF_flat.reshape((ny, nx))
        if field_component == 'efieldx':
            row = EF_2D[slice_index, :]  # Slice along x at fixed y
        else:  # efieldy
            row = EF_2D[:, slice_index]  # Slice along y at fixed x
        field_data.append(row)
    return np.vstack(field_data), len(row), time_steps

# Load data for both efieldx and efieldy
EFx, NCx, time_steps = load_field_data('efieldx', slice_index)
EFy, NCy, _ = load_field_data('efieldy', slice_index)

# Grid spacing
x = np.linspace(0, nx, nx, endpoint=False)
y = np.linspace(0, ny, ny, endpoint=False)
dx = x[1] - x[0] if nx > 1 else 1.0
dy = y[1] - y[0] if ny > 1 else 1.0

time_full = np.arange(EFx.shape[0]) * DT_coeff * write_interval

# ----------------------------- Spatial FFT and Peak Detection -----------------------------
def analyze_field(EF, NC, dx, field_name):
    # Spatial FFT
    F_k = np.fft.fft(EF, axis=1, norm='ortho')
    kfft = 2 * np.pi * np.fft.fftfreq(NC, dx)
    k_modes = kfft[:NC // 2]
    F_k_pos = F_k[:, :NC // 2]  # Positive k modes only

    # Power spectrum and peak detection
    power_spectrum = np.mean(np.abs(F_k_pos)**2, axis=0)
    peaks, _ = find_peaks(power_spectrum, height=np.max(power_spectrum) * 0.04, distance=2)

    # Harmonic analysis
    k_peaks = k_modes[peaks]
    k_fundamental = k_peaks[0] if len(peaks) > 0 else None
    harmonic_ratios = k_peaks / k_fundamental if k_fundamental is not None else []

    print(f"\n{field_name} Analysis:")
    print("Detected harmonic peak indices:", peaks)
    print("Corresponding k values:", k_modes[peaks])
    print("(k_peak / k_fundamental):", harmonic_ratios)
    print("Approximate harmonic numbers:", np.round(harmonic_ratios).astype(int))

    return F_k_pos, k_modes, power_spectrum, peaks

# Analyze both field components
Fx_k_pos, kx_modes, power_spectrum_x, peaks_x = analyze_field(EFx, NCx, dx, "Efieldx")
Fy_k_pos, ky_modes, power_spectrum_y, peaks_y = analyze_field(EFy, NCy, dy, "Efieldy")

# ----------------------------- Plot Power Spectrum with Peak Markers -----------------------------
def plot_power_spectrum(k_modes, power_spectrum, peaks, field_name, max_idx, k_peaks, harmonic_ratios):
    fig_ps, ax_ps = plt.subplots(figsize=(8, 5))
    ax_ps.plot(k_modes, power_spectrum, label='Power Spectrum')
    ax_ps.plot(k_modes[peaks], power_spectrum[peaks], 'ro', label='Detected modes')
    ax_ps.plot(k_modes[max_idx], power_spectrum[max_idx], 'ko', markersize=8,
               label=f'Max Peak at k ≈ {k_modes[max_idx]:.2f}')
    
    # Annotate fundamental and first few harmonics
    for i, (idx, ratio) in enumerate(zip(peaks[:5], harmonic_ratios[:5])):
        ax_ps.annotate(f'n={int(round(ratio))}', (k_modes[idx], power_spectrum[idx]),
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax_ps.set_xlabel('$k$')
    ax_ps.set_ylabel('$|E(k)|^2$')
    ax_ps.set_title(f'Power Spectrum ({field_name})')
    ax_ps.grid(True)
    ax_ps.legend()
    plt.tight_layout()
    plt.savefig(pjoin(path, f'power_spectrum_peaks_{field_name.lower()}.png'), dpi=300)
    plt.close()

# Plot power spectra with annotations
max_idx_x = peaks_x[np.argmax(power_spectrum_x[peaks_x])] if len(peaks_x) > 0 else 0
max_idx_y = peaks_y[np.argmax(power_spectrum_y[peaks_y])] if len(peaks_y) > 0 else 0
kx_peaks = kx_modes[peaks_x]
harmonic_ratios_x = kx_peaks / kx_peaks[0] if len(kx_peaks) > 0 else []
ky_peaks = ky_modes[peaks_y]
harmonic_ratios_y = ky_peaks / ky_peaks[0] if len(ky_peaks) > 0 else []
plot_power_spectrum(kx_modes, power_spectrum_x, peaks_x, "Efieldx", max_idx_x, kx_peaks, harmonic_ratios_x)
plot_power_spectrum(ky_modes, power_spectrum_y, peaks_y, "Efieldy", max_idx_y, ky_peaks, harmonic_ratios_y)

# ----------------------------- Energy Evolution of Detected Harmonics -----------------------------
def plot_energy_evolution(time_full, F_k_pos, k_modes, peaks, field_name):
    energy_spectrum = np.abs(F_k_pos)**2
    fig, ax = plt.subplots(figsize=(10, 6))
    num_peaks_to_plot = min(5, len(peaks))
    for idx in peaks[:num_peaks_to_plot]:
        energy = energy_spectrum[:, idx]
        ax.plot(time_full, energy, label=f'$k \\approx {k_modes[idx]:.2f}$')
    ax.set_xlabel('Time')
    ax.set_ylabel('$|E_k|^2$')
    ax.set_title(f'Time Evolution of Harmonics ({field_name})')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(pjoin(path, f'energy_exchange_between_modes_{field_name.lower()}.png'), dpi=300)
    plt.close()

# Plot energy evolution
plot_energy_evolution(time_full, Fx_k_pos, kx_modes, peaks_x, "Efieldx")
plot_energy_evolution(time_full, Fy_k_pos, ky_modes, peaks_y, "Efieldy")

# ----------------------------- Dispersion Relation Plot -----------------------------
def plot_dispersion(EF, NC, dx, time_full, field_name):
    # 2D FFT
    F = np.fft.fftn(EF, norm='ortho')
    actual_sim_time = time_full[-1] if len(time_full) > 1 else 1.0
    omega = 2 * np.pi * np.fft.fftfreq(len(time_full), time_full[1] - time_full[0] if len(time_full) > 1 else 1.0)
    k = 2 * np.pi * np.fft.fftfreq(NC, dx)
    
    # Take positive frequencies and wavenumbers
    halflen = np.array(F.shape, dtype=int) // 2
    omega = omega[:halflen[0]]
    k = k[:halflen[1]]
    F = F[:halflen[0], :halflen[1]]
    K, Omega = np.meshgrid(k, omega, indexing='xy')
    Omega /= wpe
    K *= LDe

    # Compute logarithmic power spectrum
    Z = np.abs(F)
    Z = np.log(Z + 1e-10)  # Add small constant to avoid log(0)

    # Analytic dispersion relations
    wec = e * Bz / me
    wic = e * Bz / mi
    omega_uh = np.sqrt(wpe**2 + wec**2)
    omega_lh = ((wec * wic)**(-1) + wpi**(-2))**(-0.5)
    k_vals = k * LDe
    epe = np.sqrt(wpe**2 + 3 * (wpe * LDe)**2 * k**2) / wpe  # EPW

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 5))
    c1 = ax.contourf(K, Omega, Z, cmap='rainbow', levels=50)
    plt.colorbar(c1, ax=ax, label='$\\log(|E(k,\\omega)|^2)$')
    #ax.plot(k_vals, epe, color='k', linestyle='--', lw=1.0, label='EPW')
    #ax.plot(k_vals, np.full_like(k_vals, omega_uh / wpe), color='k', linestyle='-.', lw=1.0, label='Upper Hybrid')
    #ax.plot(k_vals, np.full_like(k_vals, omega_lh / wpe), color='k', linestyle=':', lw=1.0, label='Lower Hybrid')
    ax.set_xlabel('$k \\lambda_{De}$')
    ax.set_ylabel('$\\omega / \\omega_{pe}$')
    #ax.set_xlim([0, 2])
    #ax.set_ylim([0, 2.5])
    ax.set_title(f'Dispersion Relation ({field_name})')
    ax.legend(loc='upper right', framealpha=0.5)
    plt.tight_layout()
    plt.savefig(pjoin(path, f'dispersion_{field_name.lower()}.png'), dpi=300)
    plt.close()

# Plot dispersion relations
plot_dispersion(EFx, NCx, dx, time_full, "Efieldx")
plot_dispersion(EFy, NCy, dy, time_full, "Efieldy")

# ----------------------------- kx-ky Contour Plot -----------------------------
def plot_kx_ky_contour(field_component, time_steps, time_step_idx=None):
    # Use middle time step if not specified
    if time_step_idx is None:
        time_step_idx = len(time_steps) // 2
    
    field_group = f[f'fielddata/{field_component}']
    EF_2D = field_group[str(time_steps[time_step_idx])][:].reshape((ny, nx))
    
    # Check for valid data
    if np.all(EF_2D == 0) or np.any(np.isnan(EF_2D)) or np.any(np.isinf(EF_2D)):
        print(f"Warning: {field_component} data at time step {time_steps[time_step_idx]} is invalid (all zeros, NaN, or Inf)")
        return
    
    # 2D FFT
    F_2D = np.fft.fft2(EF_2D, norm='ortho')
    kx = 2 * np.pi * np.fft.fftfreq(nx, dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, dy)
    halflen = np.array(F_2D.shape, dtype=int) // 2
    kx = kx[:halflen[1]] * LDe
    ky = ky[:halflen[0]] * LDe
    F_2D = F_2D[:halflen[0], :halflen[1]]
    Kx, Ky = np.meshgrid(kx, ky, indexing='xy')
    
    # Logarithmic power spectrum
    Z = np.abs(F_2D)
    if np.all(Z == 0):
        print(f"Warning: FFT output for {field_component} is all zeros")
        return
    Z = np.log(Z + 1e-10)  # Add small constant to avoid log(0)
    
    # Plotting
    print(f"{field_component} kx-ky plot: Z range = [{np.min(Z):.2f}, {np.max(Z):.2f}]")
    fig, ax = plt.subplots(figsize=(8, 5))
    c1 = ax.contourf(Kx, Ky, Z, cmap='rainbow', levels=50, vmin=np.percentile(Z, 5), vmax=np.percentile(Z, 95))
    plt.colorbar(c1, ax=ax, label='$\\log(|E(k_x,k_y)|^2)$')
    ax.set_xlabel('$k_x \\lambda_{De}$')
    ax.set_ylabel('$k_y \\lambda_{De}$')
    ax.set_title(f'2D Power Spectrum ({field_component}) at t={time_steps[time_step_idx]}')
    #ax.set_xlim([0, 2])
    #ax.set_ylim([0, 2])
    plt.tight_layout()
    plt.savefig(pjoin(path, f'kx_ky_spectrum_{field_component.lower()}.png'), dpi=300)
    #plt.show()
    plt.close()

# Plot kx-ky for the middle time step
plot_kx_ky_contour('efieldx', time_steps)
plot_kx_ky_contour('efieldy', time_steps)

f.close()