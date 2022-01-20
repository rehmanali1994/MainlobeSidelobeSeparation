clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd))

% Create Covariance Components
Nelem = 96; % Number of Elements on Probe
[M, S, N] = MISTcov(Nelem); % Covariance of Mainlobe, Sidelobes, and Noise

% Illustrate Covariance Components
figure; subplot(1,3,1); imagesc(M); colorbar; axis image;
xlabel('Rx Element'); ylabel('Rx Element'); title('Mainlobe Covariance');
subplot(1,3,2); imagesc(S); colorbar; axis image;
xlabel('Rx Element'); ylabel('Rx Element'); title('Sidelobes Covariance');
subplot(1,3,3); imagesc(N); colorbar; axis image;
xlabel('Rx Element'); ylabel('Rx Element'); title('Noise Covariance');

% Now take Fourier Transforms
F = dftmtx(Nelem)/Nelem; % DFT Matrix
aM = real(fftshift(diag(F*M*F'))); % FFT Bin Contributions from Mainlobe
aS = real(fftshift(diag(F*S*F'))); % FFT Bin Contributions from Sidelobes
aN = real(fftshift(diag(F*N*F'))); % FFT Bin Contributions from Noise
fftBin = fftshift(0:Nelem-1); % FFT Bin
fftBin(fftBin>=Nelem/2) = fftBin(fftBin>=Nelem/2)-Nelem; % FFT Bin Centered

% Illustrate FFT Components
figure; 
subplot(2,2,1); stem(fftBin,aM+aS,'k','Linewidth',2); 
title('Complete Aperture Spectrum');
xlabel('Rx Aperture FFT Bin'); xlim(10*[-1,1]);
ylabel('GCF Contribution'); 
subplot(2,2,2); stem(fftBin,aM,'r','Linewidth',2); 
title('Mainlobe Aperture Spectrum');
xlabel('Rx Aperture FFT Bin'); xlim(10*[-1,1]);
ylabel('GCF Contribution'); 
subplot(2,2,3); stem(fftBin,aS,'b','Linewidth',2); 
title('Sidelobe Aperture Spectrum');
xlabel('Rx Aperture FFT Bin'); xlim(10*[-1,1]);
ylabel('GCF Contribution'); 
subplot(2,2,4); stem(fftBin,aN,'g','Linewidth',2); 
title('Noise Aperture Spectrum');
xlabel('Rx Aperture FFT Bin'); xlim(10*[-1,1]);
ylabel('GCF Contribution'); 
