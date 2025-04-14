import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import os
from os.path import join as pjoin
import sys
import h5py
import matplotlib.animation as animation
import seaborn as sns
import pandas as pd

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
ne0 = n0 / (1 + alp + beta)
ni0 = n0
nn0 = alp * ne0
nb0 = beta * ne0
LD = np.sqrt(eps0 * Te / (ne0 * e**2))

sns.set_style('whitegrid')

fig, ax = plt.subplots(figsize=(10, 6))
kde_ax = fig.add_axes([0.8, 0.1, 0.1, 0.8], sharey=ax)

def animate(i):
    j = i * write_interval_phase

    data_phase_e = f["particle_electron/%d" % j]
    dataex = data_phase_e[:, 0]
    dataevx = data_phase_e[:, 1]

    data_phase_i = f["particle_ion/%d" % j]
    dataix = data_phase_i[:, 0]
    dataivx = data_phase_i[:, 1]

    #data_phase_b = f["particle_beam/%d" % j]
    #databx = data_phase_b[:, 0]
    #databvx = data_phase_b[:, 1]

    #v_total = np.concatenate((datanvx, databvx))

    ax.clear()
    kde_ax.clear()

    #sns.scatterplot(x=dataex, y= dataevx, marker='.', ax=ax, color='b', s = 10)
    sns.scatterplot(x=dataix, y= dataivx, marker='.',ax=ax, color='r', s = 10)

    #sns.histplot(y= dataevx, ax=kde_ax, color='g', kde= True, bins=100, element='step', fill=True)
    sns.histplot(y= dataivx, ax=kde_ax, color='r', kde= True, bins=100, element='step', fill=True)

    #sns.kdeplot(y= v_ion_total, ax=kde_ax, color='g', linewidth=2, fill=True)
    #sns.histplot(y=v_total, ax=kde_ax, color='g', kde= True, bins=100, element='step', fill=True)

    #sns.histplot(y=datanvx, ax=kde_ax, color='b', kde= True, bins=100, element='step', fill=True)
    #sns.histplot(y=databvx, ax=kde_ax, color='r', kde= True, bins=100, element='step', fill=True)
    #sns.violinplot(y=v_ion_total, ax=kde_ax, color='g')

    kde_ax.set_ylabel('Distribution function')
    kde_ax.get_xaxis().set_visible(False)

    return ax, kde_ax

ani = animation.FuncAnimation(fig, animate, frames=DATA_TS_PHASE, blit=False, interval=1000, repeat=False)

plt.show()
