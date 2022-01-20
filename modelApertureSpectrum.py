# Setting up all folders we can import from by adding them to python path
import sys, os, pdb

# Importing stuff from all folders in python path
import numpy as np
from scipy.ndimage import convolve
from Functions import *
import matplotlib.pyplot as plt

# Create Covariance Components
Nelem = 96; # Number of Elements on Probe
M, S, N = MISTcov(Nelem); # Covariance of Mainlobe, Sidelobes, and Noise

# Illustrate Covariance Components
plt.figure(); plt.subplot(1,3,1); imagesc(np.arange(Nelem), np.arange(Nelem), M, [np.min(M), np.max(M)]); 
plt.colorbar(fraction=0.046, pad=0.04); plt.xlabel('Rx Element'); plt.ylabel('Rx Element'); plt.title('Mainlobe Covariance');
plt.subplot(1,3,2); imagesc(np.arange(Nelem), np.arange(Nelem), S, [np.min(S), np.max(S)]); 
plt.colorbar(fraction=0.046, pad=0.04); plt.xlabel('Rx Element'); plt.ylabel('Rx Element'); plt.title('Sidelobes Covariance');
plt.subplot(1,3,3); imagesc(np.arange(Nelem), np.arange(Nelem), N, [np.min(N), np.max(N)]); 
plt.colorbar(fraction=0.046, pad=0.04); plt.xlabel('Rx Element'); plt.ylabel('Rx Element'); plt.title('Noise Covariance');
plt.tight_layout(); plt.show();

# Now take Fourier Transforms
F = np.matrix(np.fft.fft(np.eye(Nelem)))/Nelem; # DFT Matrix
aM = np.real(np.fft.fftshift(np.diag(F@M@F.H))); # FFT Bin Contributions from Mainlobe
aS = np.real(np.fft.fftshift(np.diag(F@S@F.H))); # FFT Bin Contributions from Sidelobes
aN = np.real(np.fft.fftshift(np.diag(F@N@F.H))); # FFT Bin Contributions from Noise
fftBin = np.fft.fftshift(np.arange(Nelem)); # FFT Bin
fftBin[fftBin>=Nelem/2] = fftBin[fftBin>=Nelem/2]-Nelem; # FFT Bin Centered

# Illustrate FFT Components
plt.figure(); 
plt.subplot(2,2,1); plt.stem(fftBin,aM+aS,'k',markerfmt='ko'); 
plt.title('Complete Aperture Spectrum');
plt.xlabel('Rx Aperture FFT Bin'); plt.xlim([-10,10]);
plt.ylabel('GCF Contribution'); 
plt.subplot(2,2,2); plt.stem(fftBin,aM,'r',markerfmt='ro'); 
plt.title('Mainlobe Aperture Spectrum');
plt.xlabel('Rx Aperture FFT Bin'); plt.xlim([-10,10]);
plt.ylabel('GCF Contribution'); 
plt.subplot(2,2,3); plt.stem(fftBin,aS,'b',markerfmt='bo'); 
plt.title('Sidelobe Aperture Spectrum');
plt.xlabel('Rx Aperture FFT Bin'); plt.xlim([-10,10]);
plt.ylabel('GCF Contribution'); 
plt.subplot(2,2,4); plt.stem(fftBin,aN,'g',markerfmt='go'); 
plt.title('Noise Aperture Spectrum');
plt.xlabel('Rx Aperture FFT Bin'); plt.xlim([-10,10]);
plt.ylabel('GCF Contribution'); 
plt.tight_layout(); plt.show();
