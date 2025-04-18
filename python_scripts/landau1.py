import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import value as constants
import h5py
from scipy.fft import fft, fftfreq
import sys
from os.path import join as pjoin
from scipy.optimize import curve_fit

file_name = 'result.h5'
path = sys.argv[1]

#----------------Read hdf5 file ------------
f = h5py.File(pjoin(path, file_name), 'r')
metadata_group = f['/metadata']

eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')

# Read attributes from hdf5 file 
NC = metadata_group.attrs['NC']
NUM_TS = metadata_group.attrs['NUM_TS']
write_interval  = metadata_group.attrs['write_int']
write_interval_phase = metadata_group.attrs['write_int_phase']
DT_coeff = metadata_group.attrs['DT_coeff']
nParticlesE = metadata_group.attrs['nE']
nParticlesI = metadata_group.attrs['nI']
nParticlesN = metadata_group.attrs['nN']
nParticlesB = metadata_group.attrs['nB']
Te = metadata_group.attrs['Te']
Ti = metadata_group.attrs['Ti']
Tb = metadata_group.attrs['Tb']
Tn = metadata_group.attrs['Tn']
alp = metadata_group.attrs['alpha']
beta = metadata_group.attrs['beta']
mi = metadata_group.attrs['mI']
mn = metadata_group.attrs['mN']
mb = metadata_group.attrs['mB']
n0 = metadata_group.attrs['density']
save_fig = metadata_group.attrs['save_fig']
normscheme = metadata_group.attrs['norm_scheme']
EV_TO_K = 11604.52 
DATA_TS = int(NUM_TS/write_interval) + 1

#-------------------plasma frequency and debye length calculation---------------------------
ni0 = n0
ne0 = n0/(1+alp+beta)
ni0 = n0
nn0 = alp*ne0  
nb0 = beta*ne0

LD = np.sqrt(eps0*kb*Te*EV_TO_K/(ne0*e**2)) # Characteristic Debye length
LDi = np.sqrt(eps0*kb*Ti*EV_TO_K/(ni0*e**2)) # Characteristic Debye length
we = np.sqrt(ne0*e**2/(eps0*me)) # Total Plasma Frequency
wp = np.sqrt(ni0*e**2/(eps0*mi)) # Total Plasma Frequency
DT = DT_coeff*(1.0/we)

#----------------------------------------------------------------------------------

# Transform electric field data into a 2D array
electric_field_data = []
time_steps = sorted(map(int, f['fielddata/efield'].keys()))

for i, time_step in enumerate(time_steps):
    EF_data = f['fielddata/efield/' + str(time_step)]
    electric_field_data.append(EF_data[:])  # Append electric field data

# Combine electric field data into a 2D array
EF = np.vstack(electric_field_data)

x = np.linspace(0, NC, EF.shape[1])
dx = x[1] - x[0]

wpet_1 = 0
wpet_2 = NUM_TS * DT_coeff
y1 = wpet_1 / (DT_coeff * write_interval)
y2 = wpet_2 / (DT_coeff * write_interval)

E = EF[int(y1):int(y2), :]

# Find the dominant mode amplitudes and indices
F = fft(E, axis=1, norm='ortho')
magnitude = np.abs(F)

# Find indices of top 2 maximum values
n_top_modes = 2  # Adjust for different number of dominant modes
dominant_indices = np.unravel_index(np.argsort(magnitude, axis=None)[-n_top_modes:], magnitude.shape)

# Track amplitude decay for each dominant mode
num_timesteps = E.shape[0]
amplitudes = np.zeros((num_timesteps, n_top_modes))

for t in range(num_timesteps):
    F_t = fft(E[t, :], norm='ortho')
    for k in range(n_top_modes):
        amplitudes[t, k] = np.abs(F_t[dominant_indices[k]])

# Plot amplitude decay for each dominant mode
time = np.arange(len(amplitudes)) * DT * write_interval
for k in range(n_top_modes):
    plt.plot(time, amplitudes[:, k], label=f'Dominant Mode {k+1}')

plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Amplitude Decay of the Dominant Modes')
plt.grid(True)
plt.legend()
plt.show()

# Perform exponential fit for each mode to obtain decay rates
decay_rates = []
for k in range(n_top_modes):
    popt, pcov = curve_fit(exponential_decay, time, amplitudes[:, k], p0=(amplitudes[0, k], 1.0, 0.0))
    fit_amplitudes = exponential_decay(time, *popt)
    plt.plot(time, fit_amplitudes, label=f'Exponential Fit - Mode {k+1}', linestyle='--')
    plt.legend()
    decay_rates.append(popt[1])
    print(f"Decay rate (Mode {k+1}): {decay_rates[-1]}")

plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Exponential Fit for Amplitude Decay')
plt.grid(True)
plt.show()

# Plot FFT spectrum for each dominant frequency
freqs = fftfreq(E.shape[1], d=dx)  # Frequency axis
plt.figure(figsize=(10, 6))
for k in range(n_top_modes):
    plt.plot(freqs, np.abs(F[dominant_indices[k]]), label=f'Dominant Mode {k+1}')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.title('FFT Spectrum for Dominant Frequencies')
plt.grid(True)
plt.legend()
plt.show()
