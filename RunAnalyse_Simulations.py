#####################%%%%%%
#
# Code to call simulations of a network model
# and visualize the results.
#
# Generates figure with:
# - timeseries, power spectra and phase order over time
# - matrices of coupling wights, distances and 'Functional Connectivity'
#
# Joana Cabral January 2019 joanacabral@med.uminho.pt
#
#     # Adapted to Python by Laurent Perrinet Jan 2019
#
#####################%%%%%%
import numpy as np
# Define here the Network Spatial and Temporal structure

# Define the Structural Network as a NxN matrix
from scipy.io import loadmat
C = loadmat('./AAL_matrices.mat')['C']

# Number of coupled Units
N = C.shape[0]

# note: using the line below Normalizes the whole matrix by the mean
# of all non-diagonal elements is 1. 
# C=C/mean(C(ones(N)-eye(N)>0));
# here I do what you describe in the following line:
# Normalize such that the mean of all non-diagonal elements is 1.
C[np.diag(np.ones(N))==0] /= C[np.diag(np.ones(N))==0].mean()

# Distance between areas
D = loadmat('./AAL_matrices.mat')['D']
D /= 1000 # Distance matrix in meters

# Define here the parameters of the network model to be manipulated

# Node natural frequency in Hz
f = 40    # (i.e. f=40)):

# Mean Delay in seconds
MD = 0.02 # (i.e. MD = 0.01)

# Global Coupling strength
K = 5 #


# Process the simulated data and generate figures

# Call here the function of the Network Model
from Kuramoto_Delays_Run import Kuramoto_Delays_Run_AAL
Phases_Save, dt_save = Kuramoto_Delays_Run_AAL(C, D, f, K, MD)

# Process the simulated data and generate figures
import matplotlib.pyplot as plt

N_time = Phases_Save.shape[1] 
tmax = N_time / dt_save
time = np.linspace(0, tmax, N_time)

# Plot simulated time series
fig, ax = plt.subplots(1, 1, figsize=(15, 8))
ax.plot(time, np.sin(Phases_Save.T))
ax.set_xlabel('Time (seconds)')
ax.set_ylabel('sin(\theta)')
ax.set_title(['Coupled ' + str(f) + 'Hz oscillators with for K= ' + str(K) + ' and mean delay ' + str(MD*1e3) + ' ms'])
fig.suptitle(' K= ' + str(K) + ' and mean delay ' + str(MD*1e3) + ' ms')
         
# Power Spectrum
fig, axs = plt.subplots(1, 3, figsize=(15, 8))
fbins = 5000
freqZ = np.arange(fbins)/(dt_save*fbins)
ind100Hz = np.argmin(freqZ==100)
Fourier = np.fft.fft(np.sin(Phases_Save), n=fbins, axis=-1)
PSD = np.abs(Fourier)**2
ax = axs[0]
ax.plot(freqZ[:fbins//2], PSD[:, :fbins//2].mean(axis=0))
ax.set_xlim([0, 100])
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Power')

# Plot the Order Parameter over time
ax = axs[1]
OP = np.abs(np.mean(np.exp(1j*Phases_Save), axis=0))
ax.plot(time, OP)
ax.set_ylim([0, 1])
ax.set_xlabel('Time (seconds)')
ax.set_ylabel('Order Parameter')

# Power Spectrum of the Order Parameter
ax = axs[2]
Fourier = np.fft.fft(OP-OP.mean(), n=fbins, axis=0)
PSD = np.abs(Fourier)**2
PSD = PSD[:fbins//2]
PSD /= PSD.sum()
ax.plot(freqZ[:fbins//2], PSD)
ax.set_xlim([0, 10])
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Power')

# Plot Connectivity and Distance Matrices and mean Phase Coherence

# Change the order of units for visualization (optional)
Order = np.arange(0, N, 2) # This sorts the AAL areas into left and right hemispheres
Order = np.hstack((Order, np.arange(N-1, 0, -2)))

fig, axs = plt.subplots(1, 3, figsize=(15, 8))
ax = axs[0]
# colormap(jet)
ax.pcolormesh(C[Order, :][:, Order])

ax.set_xlabel('Nodes')
ax.set_ylabel('Nodes')
ax.set_title('Coupling Matrix')
ax.axis('square')
#plt.colorbar()

ax = axs[1]
ax.pcolormesh(D[Order, :][:, Order])
ax.set_xlabel('Nodes')
ax.set_ylabel('Nodes')
ax.set_title('Distance Matrix')
ax.axis('square')
#plt.colorbar()

ax = axs[2]
X = np.sin(Phases_Save)
FC = (X[:, None, :] * X[None, :, :]).mean(axis=-1)
ax.pcolormesh(FC[Order, :][:, Order])
ax.set_xlabel('Nodes')
ax.set_ylabel('Nodes')
ax.set_title('Correlation Matrix')
ax.axis('square');
#plt.colorbar()