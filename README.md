# MainlobeSidelobeSeparation
Separation of Mainlobes and Sidelobes in the Ultrasound Image Based on the Spatial Covariance (MIST) and Aperture-Domain Spectrum of Received Signals

# Code, Results, and Sample Datasets
The underlying focusing function and covariance matrix model for MIST is implemented in both MATLAB ([Functions](Functions)) and Python ([Functions.py](Functions.py)). The following example scripts/tutorials are provided:

1) The van-Cittert Zernike theorem is used to obtain the spatial correlation between receiver signals as a function of lag (or receive element offset) in a diffuse scattering medium for the entire point spread function (PSF), the mainlobe components of the PSF, and the sidelobe components of the PSF ([theoryBehindMIST.m](theoryBehindMIST.m) and [theoryBehindMIST.py](theoryBehindMIST.py)). As suggested by the name of these scripts, these models for the spatial correlation functions of mainlobes and sidelobes in the ultrasound image lay the foundation behind multi-covariate imaging of subresolution targets (MIST). See the prior work on multicovariate imaging of subresolution targets (MIST):
> Morgan, M., Trahey, G., Walker, W. "Multi-covariate imaging of sub-resolution targets." *IEEE transactions on medical imaging* 38.7 (2019): 1690-1700.

<p align="center">
<img width="70%" height="70%" src=theoryBehindMIST.png>
</p>

2) The aforementioned spatial correlation functions are used to model the spatial covariance of received signals as well as the FFT of signals across the receive aperture ([modelApertureSpectrum.m](modelApertureSpectrum.m) and [modelApertureSpectrum.py](modelApertureSpectrum.py)). These scripts will first generate the spatial covariance matrices for the mainlobe, sidelobe, and incoherent noise contributions to the the ultrasound image as used in MIST:

<p align="center">
<img width="100%" height="100%" src=covarianceMatrices.png>
</p>

Then, these scripts will use these covariance matrices to generate the spectrum of received signals (FFT taken across the receive aperture): 

<p align="center">
<img width="90%" height="90%" src=modelApertureSpectrum.png>
</p>

3) Finally, we compare the original MIST method to our propose aperture-spectrum-based method for separating the mainlobe and sidelobe contributions to the ultrasound image ([compareReconstructedImages.m](compareReconstructedImages.m) and [compareReconstructedImages.py](compareReconstructedImages.py)). **For these specific scripts, please download the sample data (FieldII_ChannelData.mat) under the [releases](https://github.com/rehmanali1994/MainlobeSidelobeSeparation/releases) tab for this repository, and place that data in the main directory ([MainlobeSidelobeSeparation](https://github.com/rehmanali1994/MainlobeSidelobeSeparation)).**

<p align="center">
<img width="100%" height="100%" src=compareReconstructedImages.png>
</p>

# Citing this Work
If you use the code/algorithm for research, please cite the following full-length paper: 

> Rehman Ali, Trevor Mitcham, Leandra Brickson, Wentao Hu, Marvin Doyley, Deborah Rubens, Zeljko Ignjatovic, Nebojsa Duric, and Jeremy Dahl "Separation of mainlobe and sidelobe contributions to B-mode ultrasound images based on the aperture spectrum," Journal of Medical Imaging 9(6), 067001 (1 November 2022). https://doi.org/10.1117/1.JMI.9.6.067001

You can reference a static version of this code by its DOI number: [![DOI](https://zenodo.org/badge/449910520.svg)](https://zenodo.org/badge/latestdoi/449910520)
