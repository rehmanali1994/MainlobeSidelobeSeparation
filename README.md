# MainlobeSidelobeSeparation
Separation of Mainlobes and Sidelobes in the Ultrasound Image Based on the Covariance (MIST) and Aperture-Domain Spectrum of Received Signals



# Code and Sample Datasets
The underlying focusing function and covariance matrix model for MIST is implemented in both MATLAB ([Functions](Functions)) and Python ([Functions.py](Functions.py)). The following example scripts/tutorials are provided:

1) The van-Cittert Zernike theorem is used to obtain the spatial correlation between receiver signals as a function of lag (or receive element offset) in a diffuse scattering medium for the entire point spread function (PSF), the mainlobe components of the PSF, and the sidelobe components of the PSF ([theoryBehindMIST.m](theoryBehindMIST.m) and [theoryBehindMIST.py](theoryBehindMIST.py)). As suggested by the name of these scripts, these models for the spatial correlation functions of mainlobes and sidelobes in the ultrasound image lay the foundation behind multi-covariate imaging of subresolution targets (MIST). See the prior work on multicovariate imaging of subresolution targets (MIST):
> Morgan, M., Trahey, G., Walker, W. "Multi-covariate imaging of sub-resolution targets." *IEEE transactions on medical imaging* 38.7 (2019): 1690-1700.


2) The aforementioned spatial correlation functions are used to model the spatial covariance of received signals as well as the FFT of signals across the receive aperture ([modelApertureSpectrum.m](modelApertureSpectrum.m) and [modelApertureSpectrum.py](modelApertureSpectrum.py)). These scripts will first generate the spatial covariance matrices for the mainlobe, sidelobe, and incoherent noise contributions to the the ultrasound image as used in MIST:

Then, these scripts will used these covariance matrices to generate the spectrum of received signals (FFT taken across the receive aperture): 

3) Finally, we compare the original MIST method to our propose aperture-spectrum-based method for separating the mainlobe and sidelobe contributions to the ultrasound image ([compareReconstructedImages.m](compareReconstructedImages.m) and [compareReconstructedImages.py](compareReconstructedImages.py)). **Please download the sample data (FieldII_ChannelData.mat) under the [releases](https://github.com/rehmanali1994/MainlobeSidelobeSeparation/releases) tab for this repository, and place that data in the main directory ([MainlobeSidelobeSeparation](https://github.com/rehmanali1994/MainlobeSidelobeSeparation)).**

# Citing this Work
If you use the code/algorithm for research, please cite the SPIE conference proceedings paper: 

> Ali, R., Dahl, J. "Separation of Mainlobe and Sidelobe Contributions to B-Mode Ultrasound Images Based on the Aperture Spectrum". *Manuscript in preparation for SPIE Medical Imaging conference.*

You can reference a static version of this code by its DOI number: [![DOI](https://zenodo.org/badge/369032499.svg)](https://zenodo.org/badge/latestdoi/369032499)

# Results
When FWI is used to reconstruct the speed of sound in the medium using the angular spectrum method, here are the results after 12 iterations:

![](GithubResultsFigure.png)
