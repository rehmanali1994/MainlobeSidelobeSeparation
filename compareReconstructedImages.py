# Setting up all folders we can import from by adding them to python path
import sys, os, pdb

# Importing stuff from all folders in python path
import numpy as np
from scipy.ndimage import convolve
from scipy.signal import hilbert
from scipy.optimize import nnls
import scipy.io as sio
from Functions import *
import matplotlib.pyplot as plt


# Load all Functions from Subdirectories
dataset = loadmat_hdf5('FieldII_ChannelData.mat');
time = dataset['time'][0];
rxAptPos = dataset['rxAptPos'];
scat_h = hilbert(dataset['scat'], axis=0);



# Points to Focus and Get Image At
xlims = np.array([-5e-3,5e-3]); num_x = 150;
zlims = np.array([25e-3,35e-3]); num_z = 150;
x_img = np.linspace(xlims[0], xlims[1], num_x);
z_img = np.linspace(zlims[0], zlims[1], num_z);
dBrange = np.array([-60, 0]); c = 1540;

# Full Synthetic Aperture Focusing
Z, Y, X = np.meshgrid(z_img, 0, x_img);
foc_pts = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T; 
txAptPos = rxAptPos; no_elements = txAptPos.shape[0];
txApod = np.ones((num_x*num_z, no_elements));
foc_data = tx_focus_fs(time, scat_h, foc_pts, rxAptPos, txAptPos, txApod, 0, 0, c); 

# Setup MIST
[M, S, N] = MISTcov(no_elements);
MISTmodelMatrix = np.vstack((M.flatten(),S.flatten(),N.flatten())).T;


# Now take Fourier Transforms
maxBins = 6; # Values used in Proceedings Paper: 6, 12, 24, 48 (All Bins)
fftBin = np.fft.fftshift(np.arange(no_elements)); # FFT Bin
fftBin[fftBin>=no_elements/2] = fftBin[fftBin>=no_elements/2]-no_elements; # FFT Bin Centered
F = np.matrix(np.fft.fft(np.eye(no_elements)))/no_elements; # DFT Matrix
aM = np.real(np.fft.fftshift(np.diag(F@M@F.H))); # FFT Bin Contributions from Mainlobe
aS = np.real(np.fft.fftshift(np.diag(F@S@F.H))); # FFT Bin Contributions from Sidelobes
aN = np.real(np.fft.fftshift(np.diag(F@N@F.H))); # FFT Bin Contributions from Noise
fftBinModelMatrix = np.vstack((aM[abs(fftBin)<maxBins], 
    aS[abs(fftBin)<maxBins], aN[abs(fftBin)<maxBins])).T; # Model Matrix

# Focused Transmit Image Reconstruction
img_mainlobe_mist = np.zeros((num_x*num_z, 1));
img_sidelobe_mist = np.zeros((num_x*num_z, 1));
img_noise_mist = np.zeros((num_x*num_z, 1));
img_mainlobe_fft = np.zeros((num_x*num_z, 1));
img_sidelobe_fft = np.zeros((num_x*num_z, 1));
img_noise_fft = np.zeros((num_x*num_z, 1));
img_h = np.sum(foc_data,axis=1);
print('Beginning MIST and Aperture Spectrum Method');
# Create Images
for idx in np.arange(num_x*num_z):
    # Apply MIST
    obs_cov = np.outer(np.conj(foc_data[idx,:]), foc_data[idx,:]);
    asq, _ = nnls(MISTmodelMatrix, np.real(obs_cov.flatten()));
    img_mainlobe_mist[idx] = np.sqrt(asq[0]);
    img_sidelobe_mist[idx] = np.sqrt(asq[1]);
    img_noise_mist[idx] = np.sqrt(asq[2]);
    # Separate Mainlobe and Sidelobes Using Aperture Spectrum
    fftBinEnergy = np.fft.fftshift(np.abs(np.fft.fft(foc_data[idx,:]))**2).T;
    asq_fft, _ = nnls(fftBinModelMatrix, fftBinEnergy[np.abs(fftBin)<maxBins]);
    img_mainlobe_fft[idx] = np.sqrt(asq_fft[0]);
    img_sidelobe_fft[idx] = np.sqrt(asq_fft[1]);
    img_noise_fft[idx] = np.sqrt(asq_fft[2]);
img_h = img_h.reshape((num_z, num_x));
epsilon = np.max(np.abs(img_h))*np.finfo(float).eps
img_mainlobe_mist = img_mainlobe_mist.reshape((num_z, num_x))+epsilon;
img_sidelobe_mist = img_sidelobe_mist.reshape((num_z, num_x))+epsilon;
img_noise_mist = img_noise_mist.reshape((num_z, num_x))+epsilon;
img_mainlobe_fft = img_mainlobe_fft.reshape((num_z, num_x))+epsilon;
img_sidelobe_fft = img_sidelobe_fft.reshape((num_z, num_x))+epsilon;
img_noise_fft = img_noise_fft.reshape((num_z, num_x))+epsilon;

# Imaging Functions
plt.figure(); imagesc(1000*x_img, 1000*z_img, 
    20*np.log10(np.abs(img_h)/np.max(np.abs(img_h))), dBrange); 
plt.xlabel('Lateral [mm]'); plt.ylabel('Axial [mm]');
plt.title('DAS Beamforming'); plt.colorbar(); plt.show();
plt.figure(); imagesc(1000*x_img, 1000*z_img,
    20*np.log10(np.abs(img_mainlobe_mist)/np.max(np.abs(img_mainlobe_mist))), dBrange); 
plt.xlabel('Lateral [mm]'); plt.ylabel('Axial [mm]');
plt.title('MIST Mainlobe'); plt.colorbar(); plt.show();
plt.figure; imagesc(1000*x_img, 1000*z_img, 
    20*np.log10(np.abs(img_mainlobe_fft)/np.max(np.abs(img_mainlobe_fft))), dBrange); 
plt.xlabel('Lateral [mm]'); plt.ylabel('Axial [mm]');
plt.title('Aperture Spectrum Mainlobe'); plt.colorbar(); plt.show();