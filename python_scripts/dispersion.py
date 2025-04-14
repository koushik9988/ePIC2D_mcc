import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

#=================== Read Command Line Path ====================
file_name = 'result.h5'
path = sys.argv[1]

#---------------- Read hdf5 file -------------------------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

#---------------- Constants ------------------------------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#---------------- Read metadata attributes ---------------------
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']

DATA_TS = int(NUM_TS / write_interval) + 1
DT = DT_coeff * (1.0 / we)

#=================== Read 2D E-field data ======================
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efieldx'].keys()))

for time_step in time_steps:
    EF_data = f[f'fielddata/efieldx/{time_step}'][:]  # shape: (nx, ny)
    electric_field_data.append(EF_data)

EF = np.array(electric_field_data)  # shape: (time, nx, ny)
print("EF shape (t, x, y):", EF.shape)

#---------------- Space & Time Setup ---------------------------
num_ts, nx, ny = EF.shape
x = np.linspace(0, NC, nx)
y = np.linspace(0, NC, ny)
dx = x[1] - x[0]
dy = y[1] - y[0]

#================== 2D Spatial FFT =============================
EF_k = np.fft.fft2(EF, axes=(1, 2), norm='ortho')    # (t, kx, ky)
EF_k = np.fft.fftshift(EF_k, axes=(1, 2))            # shift to center

#---------------- Slice at ky = 0 ------------------------------
ky_index = ny // 2  # center corresponds to ky = 0 after fftshift
EF_kx_t = EF_k[:, :, ky_index]  # shape: (time, kx)

#================== Temporal FFT ===============================
EF_kxw = np.fft.fft(EF_kx_t, axis=0, norm='ortho')   # (omega, kx)
EF_kxw = np.fft.fftshift(EF_kxw, axes=0)             # center omega

#---------------- Frequency and Wavenumber Axes ----------------
omega = np.fft.fftshift(np.fft.fftfreq(num_ts, d=DT)) * 2 * np.pi / we
kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2 * np.pi * LD

Omega, KX = np.meshgrid(omega, kx, indexing='ij')

#================== Normalize Data for Plot ====================
Z = np.log(np.abs(EF_kxw))
#Z = np.log(np.abs(EF_kxw) / np.max(np.abs(EF_kxw)) + 1e-12)

#=================== Plot Setup ================================
figsize = np.array([80, 80 / 1.618])  # mm
dpi = 1200                            # high-res export
ppi = np.sqrt(1920**2 + 1200**2) / 24  # screen display dpi

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)

#=================== Plot Dispersion ===========================
fig, ax = plt.subplots(figsize=figsize/25.4, dpi=ppi)
c1 = ax.contourf(KX, Omega, Z, cmap='rainbow', levels=100)

cbar = plt.colorbar(c1, ax=ax)
cbar.set_label('$\zeta = \log|E(k_x, \\omega)|$')

ax.set_xlabel('$k_x \\lambda_D$')
ax.set_ylabel('$\\omega/\\omega_{pe}$')
#ax.set_xlim([0, 3])
#ax.set_ylim([0, 10])
ax.set_title("Dispersion relation (ky = 0 slice)")

#================== Save & Show ================================
plt.savefig(pjoin(path, 'dispersion_2d.png'), dpi=dpi)
plt.show()
