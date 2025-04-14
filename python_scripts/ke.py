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


data = f["time_var/kinetic_energy"]

ts = data[:,0]
keex = data[:,1]
keey = data[:,2]
keez = data[:,3]
keix = data[:,4]
keiy = data[:,5]
keiz = data[:,6]
kebx = data[:,7]
keby = data[:,8]
kebz = data[:,9]
pe = data[:,10]




fig, ((ax1, ax2),(ax3,ax4))= plt.subplots(2, 2, figsize=(10, 8))


ax1.plot(ts, keex, label='$KE_{ex}$')
ax1.plot(ts, keey,label='$KE_{ey}$')
ax1.plot(ts, keez, label='$KE_{ez}$')
ax1.set_xlabel('$\omega_{pe}t$')
ax1.set_ylabel('$KE_{e}$')
ax1.grid(True)
ax1.legend(loc='upper right',framealpha=0.5)

ax2.plot(ts, keix, label="$KE_{ix}$")
ax2.plot(ts, keiy, label="$KE_{iy}$")
ax2.plot(ts, keiz, label="$KE_{iz}$")
ax2.set_xlabel('$\omega_{pe}t$')
ax2.set_ylabel('$KE_{i}$')
ax2.grid(True)
ax2.legend(loc='upper right',framealpha=0.5)

ax3.plot(ts, kebx, label="$KE_{bx}$")
ax3.plot(ts, keby, label="$KE_{by}$")
ax3.plot(ts, kebz, label="$KE_{bz}$")
ax3.set_xlabel('$\omega_{pe}t$')
ax3.set_ylabel('$KE_{b}$')
ax3.grid(True)
ax3.legend(loc='upper right',framealpha=0.5)

ax4.plot(ts, pe, label="Total Energy")
ax4.set_xlabel('$\omega_{pe}t$')
ax4.set_ylabel('$potential energy$')
ax4.grid(True)
#ax4.set_ylim([min(pe+kee+kei) - 1.0, max(pe+kee+kei) + 1.0])
ax4.legend(loc='upper right',framealpha=0.5)

plt.tight_layout()

plt.show()