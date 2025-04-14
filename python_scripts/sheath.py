import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import time
import configparser
import h5py
import matplotlib.animation as animation


# hdf5 file name and path
file_name = 'result.h5'
path = sys.argv[1]

plot_path = './plots'

path_fig = pjoin(path, plot_path)

if not os.path.exists(path_fig):
    os.makedirs(path_fig)

# Read hdf5 file
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from input.ini file
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = e * metadata_group.attrs['Te']
Ti = e * metadata_group.attrs['Ti']
Tb = e * metadata_group.attrs['Tb']
Tn = e * metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']

DATA_TS_PHASE = int(NUM_TS / write_interval_phase) + 1

# Debye length calculation
ni0 = n0
ne0 = n0*(1-alp)
ni0 = n0
nn0 = alp * ni0

LD = np.sqrt(eps0 * Te / (ne0 * e ** 2))

j = NUM_TS
# Plotting
fig, (ax, axden) = plt.subplots(2, 1)

pot = f[f"fielddata/pot/{j}"]
den_e = f[f"fielddata/den_electron/{j}"]
den_i = f[f"fielddata/den_ion/{j}"]
den_n = f[f"fielddata/den_negion/{j}"]


x = np.linspace(0, NC, len(pot))
ax.plot(x, pot[:], color='black', label="Potential")
axden.plot(x, den_e[:], color='red', label=" electron density")
axden.plot(x, den_i[:], color='green', label="ion density")
axden.plot(x, den_n[:], color='blue', label="negative ion density")
ax.legend()
axden.legend()
ax.set_ylabel('$\phi$')
ax.legend(loc='upper right', framealpha=0.5)


plt.show()
