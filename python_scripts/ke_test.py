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
import scipy.integrate as intg


# hdf5 file name and path 
file_name = 'result.h5'
path = sys.argv[1]

path1 = './plots'

path_fig = pjoin(path, path1)

if not os.path.exists(path_fig):
  os.makedirs(path_fig)

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

# Constants and data loading from hdf5 file 
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read individual attributes
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nnParticlesN = metadata_group.attrs['nN']
nnParticlesB = metadata_group.attrs['nB']
Te = e*metadata_group.attrs['Te']
Ti = e*metadata_group.attrs['Ti']
Tb = e*metadata_group.attrs['Tb']
Tn = e*metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
EV_TO_K = 11604.52 

# Get number of particle species from user input
n_species = int(input("Enter the number of species: "))

# Define empty lists to store particle data
ke_data = []
species_names = []

# Loop through each species and load data
for i in range(n_species):
  species_name = f'species{i+1}'  # Assuming species names are like 'species1', 'species2' etc.
  species_names.append(species_name)
  ke_data.append(f["time_var/kinetic_energy_" + species_name][:,1:])  # Assuming kinetic energy data starts from second column

# Debye lenght calculation (assuming single electron species)
ni0 = n0
ne0 = n0*((1-alp-beta*beta))
ni0 = n0
nn0 = alp*ni0  
nb0 = beta*ni0

LD = np.sqrt(eps0*Te/(ne0*e**2)) 
electron_spwt = (ne0*NC*LD)/(nParticlesE)
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency

#---------------------------------------------------------
data = f["time_var/kinetic_energy"]

ts = data[:,0]


#-----potential-energy calculation----------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
PE = np.zeros(len(time_steps))
#print(time_steps)
for i, time_step in enumerate(time_steps):
  EF_data = f['fielddata/efield/' + str(time_step)]
  x = np.linspace(0,NC,len(EF_data))*LD
  E_sq = (EF_data[:]** 2) * ((me*we**2*LD/e)**2) 
  integral = intg.trapz(E_sq, x)
  PE[i] = 0.5 * eps0 * integral

THe = ((electron_spwt)*nParticlesE)*Te#*EV_TO_K*kb
PE/= THe
#---------------plotting -------------------------

num_rows = (n_species + 2) // 2  # Calculate number of rows for subplot
num_cols = min(n_species, 3)  # Limit columns to maximum 3

fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))  # Adjust figsize as needed

# Plotting kinetic energy for each species
for i, ke in enumerate(ke_data):
  row = i // num_cols
  col = i % num_cols
  axes[row, col].plot(ts, ke, label=species_names[i])
