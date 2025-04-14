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


#hdf5 file name and path 
file_name = 'result.h5'
path = sys.argv[1]

path1 = './plots'

path_fig = pjoin(path,path1)

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
normscheme = metadata_group.attrs['norm_scheme']
EV_TO_K = 11604.52 

#-----------------------Debye lenght calculation--------
ni0 = n0
ne0 = n0/(1+alp+beta)
ni0 = n0
nn0 = alp*ne0  
nb0 = beta*ne0

# Characteristic Debye length
LD = np.sqrt(eps0*Te/(ne0*e**2)) 
LDi = np.sqrt(eps0*Ti/(ni0*e**2)) 
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency
#---------------------------------------------------------
data = f["time_var/kinetic_energy"]

ts = data[:,0]
kee = data[:,1]
kei = data[:,2]
ken = data[:,3]
keb = data[:,4]


if normscheme == 2 or normscheme == 4:
    LD = LDi
if normscheme == 1 or normscheme == 3:
    LD = LD

electron_spwt = (ne0*NC*LD)/(nParticlesE)
ion_spwt = (ni0*NC*LD)/(nParticlesI)
print(electron_spwt)
#-----potential-energy calculation----------
time_steps = sorted(map(int, f['fielddata/efield'].keys()))
PE = np.zeros(len(time_steps))
#print(time_steps)
for i, time_step in enumerate(time_steps):
        EF_data = f['fielddata/efield/' + str(time_step)]
        x = np.linspace(0,NC,len(EF_data))*LD
        #E_sq = (EF_data[:]** 2) * ((me*we**2*LD/e)**2)
        E_sq = (EF_data[:]** 2) * ((ni0*e*LD)/eps0)**2  
        integral = intg.trapz(E_sq, x)
        PE[i] = 0.5 * eps0 * integral

THe = ((electron_spwt)*nParticlesE)*Te#*EV_TO_K*kb

THi = ((ion_spwt)*nParticlesI)*Ti#*EV_TO_K*kb

#normalize according to normalization 
if normscheme == 2 or normscheme == 4:
     PE/= THi
if normscheme == 1 or normscheme == 3:
    PE /= THe
#---------------plotting -------------------------

fig, ((ax1, ax2, ax5),(ax6, ax3, ax4))= plt.subplots(2, 3, figsize=(10, 8))


ax1.plot(ts, kee, label='$KE_{e}$')
ax1.set_xlabel('$\omega_{pe}t$')
ax1.set_ylabel('$KE_{e}$')
ax1.grid(True)
ax1.legend(loc='upper right',framealpha=0.5)

ax2.plot(ts, kei, label="$KE_{i}$")
ax2.set_xlabel('$\omega_{pe}t$')
ax2.set_ylabel('$KE_{i}$')
ax2.grid(True)
ax2.legend(loc='upper right',framealpha=0.5)

ax5.plot(ts, ken, label="$KE_{n}$")
ax5.set_xlabel('$\omega_{pe}t$')
ax5.set_ylabel('$KE_{n}$')
ax5.grid(True)
ax5.legend(loc='upper right',framealpha=0.5)

ax6.plot(ts, keb , label="$KE_{b}$")
ax6.set_xlabel('$\omega_{pe}t$')
ax6.set_ylabel('$KE_{b}$')
ax6.grid(True)
ax6.legend(loc='upper right',framealpha=0.5)

ax3.plot(ts, PE, label="potential energy")
ax3.set_xlabel('$\omega_{pe}t$')
ax3.set_ylabel('$PE$')
ax3.grid(True)
ax3.legend(loc='upper right',framealpha=0.5)

ax4.plot(ts, kei + kee + keb + ken + PE, label="Total Energy")
ax4.set_xlabel('$\omega_{pe}t$')
ax4.set_ylabel('Total Energy')
ax4.grid(True)
ax4.set_ylim([min(kei + kee + keb + ken + PE) - 1.0, max(kei + kee + keb + ken + PE) + 1.0])
ax4.legend(loc='upper right',framealpha=0.5)

plt.tight_layout()

#if(save_fig == 1):
plt.savefig(pjoin(path,'ke_pe.png'),dpi = 1200)

plt.show()