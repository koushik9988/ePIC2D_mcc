import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import h5py
import sys
from os.path import join as pjoin

#------ Load Arguments and Constants ------
path = sys.argv[1]
file_name = 'result.h5'

eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#------ Open HDF5 Data ------
f = h5py.File(pjoin(path, file_name), 'r')
meta = f['/metadata']

NC = meta.attrs['NC']
NUM_TS = meta.attrs['NUM_TS']
DT_coeff = meta.attrs['DT_coeff']
write_interval = meta.attrs['write_int']
LD = meta.attrs['LDe']
we = meta.attrs['wpe']

DT = DT_coeff / we  # Time step in normalized units
DATA_TS = int(NUM_TS / write_interval) + 1

#------ Load Electric Field Data ------
efield_data = []
time_steps = sorted(map(int, f['fielddata/efieldx'].keys()))

for step in time_steps:
    E = f[f'fielddata/efieldy/{step}'][:].reshape((NC, NC))
    efield_data.append(E)

EF = np.stack(efield_data, axis=0)  # Shape: (time, nx, ny)
print(f"Loaded E-field shape: {EF.shape}")

#------ Perform FFT ------
F = np.fft.fftn(EF, axes=(0, 1, 2), norm='ortho')
F = np.fft.fftshift(F, axes=(0, 1, 2))

#------ Frequency & Wavenumber Axes ------
dt = DT
dx = 1.0
dy = 1.0

omega = np.fft.fftshift(np.fft.fftfreq(EF.shape[0], d=dt)) * 2 * np.pi / we
kx = np.fft.fftshift(np.fft.fftfreq(EF.shape[1], d=dx)) * 2 * np.pi * LD
ky = np.fft.fftshift(np.fft.fftfreq(EF.shape[2], d=dy)) * 2 * np.pi * LD

#------ Select ky Slice ------
# Choose ky=0 mode first:
ky_target = 0
ky_idx = np.abs(ky - ky_target).argmin()

# Extract and log-scale
Z_slice = np.log(np.abs(F[:, :, ky_idx]) + 1e-20)

#------ Plotting ------
fig, ax = plt.subplots()
omega_mesh, kx_mesh = np.meshgrid(omega, kx, indexing='ij')

contour = ax.contourf(kx_mesh, omega_mesh, Z_slice, cmap='rainbow', levels=100)
plt.colorbar(contour, label=r'$\log|\tilde{E}(k_x, \omega)|$')

ax.set_xlabel(r'$k_x \lambda_D$')
ax.set_ylabel(r'$\omega / \omega_{pe}$')
ax.set_title(f"Dispersion Plot at $k_y = {ky[ky_idx]:.2f}$")

plt.savefig(pjoin(path, 'dispersion_2D.png'), dpi=300)
plt.show()
