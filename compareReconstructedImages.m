clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd))

% Load Channel Data
load('FieldII_ChannelData.mat');

% Points to Focus and Get Image At
xlims = [-5e-3,5e-3]; num_x = 150;
zlims = [25e-3,35e-3]; num_z = 150;
x_img = linspace(xlims(1), xlims(2), num_x);
z_img = linspace(zlims(1), zlims(2), num_z);
dBrange = [-60, 0]; c = 1540;

% Full Synthetic Aperture Focusing
[Z, Y, X] = meshgrid(z_img, 0, x_img);
foc_pts = [X(:), Y(:), Z(:)]; 
txAptPos = rxAptPos; no_elements = size(txAptPos,1);
txApod = ones(num_x*num_z, no_elements);
scat_h = hilbert(scat);
tic; foc_data = tx_focus_fs_fast(time, scat_h, ...
    foc_pts, rxAptPos, txAptPos, txApod, 0, 0, c); toc;

% Setup MIST
[M, S, N] = MISTcov(no_elements);
MISTmodelMatrix = [M(:),S(:),N(:)];

% Now take Fourier Transforms
maxBins = 6; % Values used in Proceedings Paper: 6, 12, 24, 48 (All Bins)
fftBin = fftshift(0:no_elements-1); % FFT Bin
fftBin(fftBin>=no_elements/2) = ...
    fftBin(fftBin>=no_elements/2)-no_elements; % FFT Bin Centered
F = dftmtx(no_elements)/no_elements; % DFT Matrix
aM = real(fftshift(diag(F*M*F'))); % FFT Bin Contributions from Mainlobe
aS = real(fftshift(diag(F*S*F'))); % FFT Bin Contributions from Sidelobes
aN = real(fftshift(diag(F*N*F'))); % FFT Bin Contributions from Noise
fftBinModelMatrix = [aM(abs(fftBin)<maxBins), ...
    aS(abs(fftBin)<maxBins), aN(abs(fftBin)<maxBins)]; % Model Matrix

% Focused Transmit Image Reconstruction
img_mainlobe_mist = zeros(num_x*num_z, 1);
img_sidelobe_mist = zeros(num_x*num_z, 1);
img_noise_mist = zeros(num_x*num_z, 1);
img_mainlobe_fft = zeros(num_x*num_z, 1);
img_sidelobe_fft = zeros(num_x*num_z, 1);
img_noise_fft = zeros(num_x*num_z, 1);
img_h = sum(foc_data,2);
disp('Beginning MIST and Aperture Spectrum Method');
% Create Images
for idx = 1:num_x*num_z
    % Apply MIST
    obs_cov = foc_data(idx,:)'*foc_data(idx,:);
    asq = lsqnonneg(MISTmodelMatrix, real(obs_cov(:)));
    img_mainlobe_mist(idx) = sqrt(asq(1));
    img_sidelobe_mist(idx) = sqrt(asq(2));
    img_noise_mist(idx) = sqrt(asq(3));
    % Separate Mainlobe and Sidelobes Using Aperture Spectrum
    fftBinEnergy = fftshift(abs(fft(foc_data(idx,:))).^2)';
    asq_fft = lsqnonneg(fftBinModelMatrix, fftBinEnergy(abs(fftBin)<maxBins));
    img_mainlobe_fft(idx) = sqrt(asq_fft(1));
    img_sidelobe_fft(idx) = sqrt(asq_fft(2));
    img_noise_fft(idx) = sqrt(asq_fft(3));
end
img_h = reshape(img_h, [num_z, num_x]);
img_mainlobe_mist = reshape(img_mainlobe_mist, [num_z, num_x]);
img_sidelobe_mist = reshape(img_sidelobe_mist, [num_z, num_x]);
img_noise_mist = reshape(img_noise_mist, [num_z, num_x]);
img_mainlobe_fft = reshape(img_mainlobe_fft, [num_z, num_x]);
img_sidelobe_fft = reshape(img_sidelobe_fft, [num_z, num_x]);
img_noise_fft = reshape(img_noise_fft, [num_z, num_x]);

% Imaging Functions
figure; imagesc(1000*x_img, 1000*z_img, ...
    20*log10(abs(img_h)/max(abs(img_h(:)))), dBrange); 
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('DAS Beamforming'); colormap(gray); colorbar();
figure; imagesc(1000*x_img, 1000*z_img, ...
    20*log10(abs(img_mainlobe_mist)/max(abs(img_mainlobe_mist(:)))), dBrange); 
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('MIST Mainlobe'); colormap(gray); colorbar();
figure; imagesc(1000*x_img, 1000*z_img, ...
    20*log10(abs(img_mainlobe_fft)/max(abs(img_mainlobe_fft(:)))), dBrange); 
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('Aperture Spectrum Mainlobe'); colormap(gray); colorbar();