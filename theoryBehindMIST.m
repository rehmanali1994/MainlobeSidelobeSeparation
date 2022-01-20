clear
clc

% Parameters
Nelem = 96; % Number of Elements on Probe
Fs = 10000; % Sampling Rate (Samples Per Element)--Large Number for Pseudoanalytic Solution

% Spatial and Wavenumber Axes
Nx = Fs*Nelem+1;
x = linspace(-Nelem/2,Nelem/2,Nx); % Aperture Coordinate (in Elements)
kx = Fs*(-(Nx-1)/2:(Nx-1)/2)/Nx; % Arbitrary Lateral Coordinate for Far-Field Pattern

% Create Wavenumber-Domain Signals and Their Inverse Fourier Transforms
psf_vcz_orig = max(0,1-abs(kx)); % Linearly Decreasing VCZ curve
box_vcz_orig = 2*sinc(2*kx); % Fourier Transform of Rect Function over the Mainlobe
mainlobe_vcz_orig = conv(psf_vcz_orig,box_vcz_orig)/Nelem; % Convolve to Get Mainlobe Only
mainlobe_vcz_orig = mainlobe_vcz_orig((-(Nx-1)/2:(Nx-1)/2)+Nx);
psf = fftshift(fft(ifftshift(psf_vcz_orig)))/Nelem; % Full PSF
box = fftshift(fft(ifftshift(box_vcz_orig)))/Nelem; % Box Over the Mainlobe
mainlobe = fftshift(fft(ifftshift(mainlobe_vcz_orig)))/Nelem; % Mainlobe of the PSF
sidelobe = psf-mainlobe; % Sidelobes of the PSF

% Generate VCZ Curves for Mainlobe and Sidelobes 
xcorr_aperture = (kx>=0 & kx<=1);
lag = kx(xcorr_aperture);
full_vcz = real(psf_vcz_orig(xcorr_aperture));
mnlobe_vcz = real(mainlobe_vcz_orig(xcorr_aperture));
sdlobe_vcz = full_vcz-mnlobe_vcz;
figure; subplot(2,1,1); plot(kx, psf, ...
    kx, mainlobe, '--', kx, sidelobe, '--', 'Linewidth', 2);
legend('Full PSF', 'Mainlobe Only', 'Sidelobe Only');
xlabel('Lateral Coordinate'); ylabel('Intensity'); xlim(3*[-Nelem,Nelem]); 
set(gca,'xticklabel',{[]}); set(gca,'yticklabel',{[]});
title('Components of Point Spread Function (PSF)');
subplot(2,1,2); plot(lag, full_vcz, ...
    lag, mnlobe_vcz, lag, sdlobe_vcz, 'Linewidth', 2);
legend('Full VCZ', 'Mainlobe Only', 'Sidelobe Only');
xlabel('lag (fraction of aperture)'); ylabel('Correlation');
title('Van-Cittert Zernike (VCZ) Theorem for PSF Components');