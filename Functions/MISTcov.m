function [mainlobe_cov, sidelobe_cov, noise_cov] = MISTcov(nRx)
%MISTcov Generate Covariance Matrices for Mainlobe, Sidelobes, and Noise
%   Inputs:
%       nRx = number of recieve channels
%   Outputs:
%       mainlobe_cov = covariance matrix for mainlobe
%       sidelobe_cov = covariance matrix for sidelobes
%       noise_cov = covariance matrix for noise

% Sampling Rate (Samples Per Element)
Fs = 10000; % Make this a Large Number for Pseudoanalytic Solution

% Spatial and Wavenumber Axes
Nx = Fs*nRx+1; % Total Number of X Samples
kx = Fs*(-(Nx-1)/2:(Nx-1)/2)/Nx; % Arbitrary Lateral Coordinate for Far-Field Pattern

% Create Wavenumber-Domain Signals and Their Inverse Fourier Transforms
psf_vcz_orig = max(0,1-abs(kx)); % Linearly Decreasing VCZ curve
box_vcz_orig = 2*sinc(2*kx); % Fourier Transform of Rect Function over the Mainlobe
mainlobe_vcz_orig = conv(psf_vcz_orig,box_vcz_orig)/nRx; % Convolve to Get Mainlobe Only
mainlobe_vcz_orig = mainlobe_vcz_orig((-(Nx-1)/2:(Nx-1)/2)+Nx);

% Generate VCZ Curves for Mainlobe and Sidelobes 
xcorr_aperture = (kx>=0 & kx<=1);
lag = kx(xcorr_aperture);
full_vcz = real(psf_vcz_orig(xcorr_aperture));
mnlobe_vcz = real(mainlobe_vcz_orig(xcorr_aperture));
sdlobe_vcz = full_vcz-mnlobe_vcz;

% Form Covariance Matrices for Mainlobe, Sidelobes, and Noise
sidelobe_vcz = interp1(lag, sdlobe_vcz, (0:nRx-1)/nRx);
sidelobe_cov = toeplitz(sidelobe_vcz, sidelobe_vcz);
mainlobe_vcz = interp1(lag, mnlobe_vcz, (0:nRx-1)/nRx);
mainlobe_cov = toeplitz(mainlobe_vcz, mainlobe_vcz);
noise_cov = eye(nRx); % Covariance Matrix for Noise

end