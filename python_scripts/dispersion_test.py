import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import h5py

# === Input Arguments ===
if len(sys.argv) < 4:
    print("Usage: python script.py <path> <field_component: efieldx|efieldy> <slice_index>")
    sys.exit(1)

path = sys.argv[1]
field_component = sys.argv[2].strip().lower()
slice_index = int(sys.argv[3])

if field_component not in ['efieldx', 'efieldy']:
    print("Error: field_component must be 'efieldx' or 'efieldy'")
    sys.exit(1)

file_name = 'result.h5'
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']
metadata_electron = f['/metadata_species/electron']
metadata_ion = f['/metadata_species/ion']

# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Metadata
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
DT_coeff = metadata_group.attrs['DT_coeff']
LD = metadata_group.attrs['LDe']
LDi = metadata_group.attrs['LDi']
we = metadata_group.attrs['wpe']
wp = metadata_group.attrs['wpi']
n0 = metadata_group.attrs['density']
normscheme = metadata_group.attrs['norm_scheme']
Te = metadata_electron.attrs['temperature']*e
Ti = metadata_ion.attrs['temperature']*e
mi = metadata_ion.attrs['mass']


# Derived Quantities
vthi = np.sqrt(Ti / mi)
vthe = np.sqrt(Te / me)

mfactor = 1 if normscheme in [1, 2] else wp / we
DT = DT_coeff * (1.0 / we)

# Grid Size
nx = metadata_group.attrs['Nx']
ny = metadata_group.attrs['Ny']
NC = nx if field_component == 'efieldx' else ny

# Load field data and extract slice
efield_group = f[f'fielddata/{field_component}']
time_steps = sorted(map(int, efield_group.keys()))
electric_field_data = []

for time_step in time_steps:
    EF_flat = efield_group[str(time_step)][:]
    EF_2D = EF_flat.reshape((ny, nx))
    
    if field_component == 'efieldx':
        if slice_index >= ny:
            print(f"Error: slice_index out of bounds for ny={ny}")
            sys.exit(1)
        row = EF_2D[slice_index, :]
    else:  # efieldy
        if slice_index >= nx:
            print(f"Error: slice_index out of bounds for nx={nx}")
            sys.exit(1)
        row = EF_2D[:, slice_index]
    
    electric_field_data.append(row)

EF = np.vstack(electric_field_data) 
dx = 1.0 

wpet_1 = 0
wpet_2 = NUM_TS * DT_coeff
y1 = int(wpet_1 / (DT_coeff * write_interval))
y2 = int(wpet_2 / (DT_coeff * write_interval))
E = EF[y1:y2, :]

# FFT
F = np.fft.fftn(E, norm='ortho')
NUM_TS2 = wpet_2 / DT_coeff
NUM_TS1 = wpet_1 / DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1) * DT * mfactor
omega = 2 * np.pi * np.arange(E.shape[0]) / actual_sim_time
k = 2 * np.pi * np.arange(NC + 1) / (NC * dx * LD)

Omega, K = np.meshgrid(omega, k, indexing='ij')
halflen = np.array(F.shape, dtype=int) // 2
Omega = Omega[:halflen[0], :halflen[1]]
K = K[:halflen[0], :halflen[1]]
F = F[:halflen[0], :halflen[1]]
Omega /= we
K *= LD

Z = np.log(np.abs(F))


Bz = 0.0001

wec = e*Bz/me

wic = e*Bz/mi

omega_uh = np.sqrt(we**2 + wec**2)

omega_lh = ((wec*wic)**(-1) + wp**(-2))**(-0.5)

# Analytic EPW
raw_analytic_EPW = True
if raw_analytic_EPW:
    epe = np.sqrt(we**2 + 3 * (we * LD)**2 * k**2)

raw_analytic_upperhybrid = True
if raw_analytic_EPW:
    wuh = np.full(len(k), omega_uh)

raw_analytic_lowerhybrid = True
if raw_analytic_EPW:
    wlh = np.full(len(k), omega_lh)

# Plotting
figsize = np.array([80, 80 / 1.618])
dpi = 300
ppi = np.sqrt(1920**2 + 1200**2) / 24
mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10)

fig, ax = plt.subplots(figsize=figsize / 10.4, constrained_layout=True, dpi=ppi)
c1 = plt.contourf(K, Omega, Z, cmap='rainbow', shading='auto')#, vmin= 0, vmax= 2.5)
#extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
#c1 = plt.imshow(Z, origin='lower', interpolation='bilinear', cmap='coolwarm')
#cbar = plt.colorbar()
#cbar.set_label('$\zeta$')

ax.set_xlabel('$k \lambda_{De}$')
ax.set_ylabel('$\omega/\omega_{pe}$')
ax.set_xlim([0, 2])
ax.set_ylim([0, 2.5])

if raw_analytic_EPW:
    plt.plot(k * LD, epe / we, color='k', linestyle='--', lw=1.0, label='EPW')

if raw_analytic_upperhybrid:
     plt.plot(k*LD, wuh/we, color='k', linestyle='--', lw = 1.0, label='$upper\_hybrid$')

if raw_analytic_lowerhybrid:
     plt.plot(k*LD, wlh/we, color='k', linestyle='-.', lw = 1.0, label='$lower\_hybrid$')

ax.legend(loc='upper right', framealpha=0.5)
plt.savefig(pjoin(path, 'dispersion.png'), dpi=dpi)
plt.show()
