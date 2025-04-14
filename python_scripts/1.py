import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import h5py
import sys
from os.path import join as pjoin

#----------------- Arguments ------------------
file_name = 'result.h5'
path = sys.argv[1]

#----------------- HDF5 File ------------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

#---------------- Constants -------------------
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

#--------------- Attributes -------------------
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
DT_coeff = metadata_group.attrs['DT_coeff']
write_interval = metadata_group.attrs['write_int']
LD = metadata_group.attrs['LDe']
we = metadata_group.attrs['wpe']
normscheme = metadata_group.attrs['norm_scheme']

#--------------- Time Step and Scaling --------
DT = DT_coeff / we
DATA_TS = int(NUM_TS / write_interval) + 1

#------------- Load E-field -------------------
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efieldx'].keys()))

for time_step in time_steps:
    EF_data = f[f'fielddata/den_electron/{time_step}']
    EF = np.array(EF_data[:])
    EF = EF.reshape((metadata_group.attrs['NC'], metadata_group.attrs['NC']))  # from 1D to 2D
    electric_field_data.append(EF)

EF = np.stack(electric_field_data, axis=0)  # Shape: (time, nx, ny)
print("EF shape (time, nx, ny):", EF.shape)

#------------ FFT -----------------------------
F = np.fft.fftn(EF, axes=(0, 1, 2), norm='ortho')  # full 3D FFT
F = np.fft.fftshift(F, axes=(0, 1, 2))             # shift zero-frequency to center

# Only positive frequency half for plotting
half_t = EF.shape[0] // 2
half_x = EF.shape[1] // 2
half_y = EF.shape[2] // 2
F = F[half_t:, half_x:, half_y:]

Z = np.log(np.abs(F) + 1e-20)  # avoid log(0)

#------ Frequency and Wavenumber Axes ---------
t_steps = EF.shape[0]
nx = EF.shape[1]
ny = EF.shape[2]
dx = 1.0  # update with actual dx if needed
dy = 1.0  # update with actual dy if needed

omega = np.fft.fftshift(np.fft.fftfreq(t_steps, d=DT)) * 2 * np.pi
kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2 * np.pi * LD
ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy)) * 2 * np.pi * LD

omega = omega[half_t:] / we
kx = kx[half_x:]
ky = ky[half_y:]

#--------- Plot Settings -----------------------
fig, ax = plt.subplots()
#X, Y = np.meshgrid(kx, omega, indexing='ij')  # or use ky if plotting vs ky
Y, X = np.meshgrid(omega, kx, indexing='ij')

Z_slice = Z[:, :, Z.shape[2] // 2]  # take central slice in ky (or kx)

c = ax.contourf(X, Y, Z_slice, cmap='rainbow', levels=100)
cbar = plt.colorbar(c)
cbar.set_label(r'$\log|\tilde{E}(k,\omega)|$')

ax.set_xlabel(r'$k_x \lambda_D$')
ax.set_ylabel(r'$\omega/\omega_{pe}$')
#ax.set_xlim([0, 3.0])
#ax.set_ylim([0, 10])
ax.set_title("Dispersion Relation (central slice in $k_y$)")

#----------- Save and Show ---------------------
plt.savefig(pjoin(path, 'dispersion_2D.png'), dpi=300)
plt.show()
