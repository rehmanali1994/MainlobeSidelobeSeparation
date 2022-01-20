# Setting up all folders we can import from by adding them to python path
import sys, os, pdb

# Importing stuff from all folders in python path
import numpy as np
from scipy.ndimage import convolve
from Functions import *
import matplotlib.pyplot as plt

# Parameters
Nelem = 96; # Number of Elements on Probe
Fs = 10000; # Sampling Rate (Samples Per Element)--Large Number for Pseudoanalytic Solution

# Spatial and Wavenumber Axes
Nx = Fs*Nelem+1;
x = np.linspace(-Nelem/2,Nelem/2,Nx); # Aperture Coordinate (in Elements)
kx = Fs*np.arange(-(Nx-1)/2,1+(Nx-1)/2)/Nx; # Arbitrary Lateral Coordinate for Far-Field Pattern

# Create Wavenumber-Domain Signals and Their Inverse Fourier Transforms
psf_vcz_orig = np.maximum(0,1-np.abs(kx)); # Linearly Decreasing VCZ curve
box_vcz_orig = 2*np.sinc(2*kx); # Fourier Transform of Rect Function over the Mainlobe
mainlobe_vcz_orig = np.real(np.fft.ifft(np.fft.fft(psf_vcz_orig, n=2*(Nx-1)+1)* 
    np.fft.fft(box_vcz_orig, n=2*(Nx-1)+1)))/Nelem; # Convolve to Get Mainlobe Only
mainlobe_vcz_orig = mainlobe_vcz_orig[int(Nx-1-(Nx-1)/2):int(Nx+(Nx-1)/2)];
psf = np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(psf_vcz_orig))))/Nelem; # Full PSF
box = np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(box_vcz_orig))))/Nelem; # Box Over the Mainlobe
mainlobe = np.real(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(mainlobe_vcz_orig))))/Nelem; # Mainlobe of the PSF
sidelobe = psf-mainlobe; # Sidelobes of the PSF

# Generate VCZ Curves for Mainlobe and Sidelobes 
xcorr_aperture = np.logical_and(kx>=0, kx<=1);
lag = kx[xcorr_aperture];
full_vcz = np.real(psf_vcz_orig[xcorr_aperture]);
mnlobe_vcz = np.real(mainlobe_vcz_orig[xcorr_aperture]);
sdlobe_vcz = full_vcz-mnlobe_vcz;
plt.figure(); plt.subplot(2,1,1); 
plt.plot(kx, psf, label = 'Full PSF', linewidth = 2)
plt.plot(kx, mainlobe, '--', label = 'Mainlobe Only', linewidth = 2);
plt.plot(kx, sidelobe, '--', label = 'Sidelobe Only', linewidth = 2);
plt.legend(); ax = plt.gca(); ax.axes.xaxis.set_ticklabels([]); ax.axes.yaxis.set_ticklabels([]);
plt.xlabel('Lateral Coordinate'); plt.ylabel('Intensity'); plt.xlim([-3*Nelem,3*Nelem]);
plt.title('Components of Point Spread Function (PSF)');
plt.subplot(2,1,2); plt.plot(lag, full_vcz, label = 'Full VCZ', linewidth = 2)
plt.plot(lag, mnlobe_vcz, label = 'Mainlobe Only', linewidth = 2);
plt.plot(lag, sdlobe_vcz, label = 'Sidelobe Only', linewidth = 2);
plt.xlabel('lag (fraction of aperture)'); plt.ylabel('Correlation');
plt.title('Van-Cittert Zernike (VCZ) Theorem for PSF Components');
plt.tight_layout(); plt.xlim([0,1]); plt.show();