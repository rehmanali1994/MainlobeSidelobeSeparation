# MainlobeSidelobeSeparation
Separation of Mainlobes and Sidelobes in the Ultrasound Image Based on the Covariance (MIST) and Aperture-Domain Spectrum of Received Signals



# Code and Sample Datasets
The underlying focusing function and covariance matrix model for MIST is implemented in both MATLAB ([Functions](Functions)) and Python ([Functions.py](Functions.py)). The following example scripts/tutorials are provided:

1) The van-Cittert Zernike theorem is used to obtain the spatial correlation between receiver signals as a function of lag (or receive element offset) in a diffuse scattering medium for the entire point spread function (PSF), the mainlobe components of the PSF, and the sidelobe components of the PSF ([theoryBehindMIST.m](theoryBehindMIST.m) and [theoryBehindMIST.py](theoryBehindMIST.py)). As suggested by the name of these scripts, these models for spatial correlation function lay the foundation behind multi-covariate imaging of subresolution targets (MIST). See the prior work on multicovariate imaging of subresolution targets (MIST):
> Morgan, M., Trahey, G., Walker, W. "Multi-covariate imaging of sub-resolution targets." *IEEE transactions on medical imaging* 38.7 (2019): 1690-1700.

2) The equivalent time-domain reconstruction process shown using a single-element transmission ([TimeDomShotGatherMig_FieldII.m](TimeDomShotGatherMig_FieldII.m) and [TimeDomShotGatherMig_FieldII.py](TimeDomShotGatherMig_FieldII.py)).
3) A focused synthetic aperture image reconstruction for an abdominal imaging example ([FreqDomShotGatherMig_Curvilinear_5C1.m](FreqDomShotGatherMig_Curvilinear_5C1.m) and [FreqDomShotGatherMig_Curvilinear_5C1.py](FreqDomShotGatherMig_Curvilinear_5C1.py)). Special thank to **Rick Loftman**, **Ismayil Guracar**, and **Vasant Salgaonkar** at **Siemens Healthineers** for enabling in-vivo demonstration of this work using channel data captured on a clinical scanner.

**Please download the sample data (FieldII_ChannelData.mat) under the [releases](https://github.com/rehmanali1994/MainlobeSidelobeSeparation/releases) tab for this repository, and place that data in the main directory ([MainlobeSidelobeSeparation](https://github.com/rehmanali1994/MainlobeSidelobeSeparation)).**

# Citing this Work
If you use the code/algorithm for research, please cite the SPIE conference proceedings paper: 

> Ali, R., Dahl, J. "Separation of Mainlobe and Sidelobe Contributions to B-Mode Ultrasound Images Based on the Aperture Spectrum". *Manuscript in preparation for SPIE Medical Imaging conference.*

You can reference a static version of this code by its DOI number: [![DOI](https://zenodo.org/badge/369032499.svg)](https://zenodo.org/badge/latestdoi/369032499)

# Schematic of the Imaging System
The schematic below shows the coordinate system (a) used to perform the angular spectrum method. This grid is rotated as the two linear arrays are (b) rotated around the medium to collect receive signals sampled from all angles by rotating a full a 360 degrees (in 2 degree steps) around the object of interest:

![](TxTomography.png)

# Simulated Dataset
Simulated receive signals (from k-Wave) are shown for 9 different views around the medium:

![](FullWaveformSetup.png)

# Results
When FWI is used to reconstruct the speed of sound in the medium using the angular spectrum method, here are the results after 12 iterations:

![](GithubResultsFigure.png)
