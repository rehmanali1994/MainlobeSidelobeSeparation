import numpy as np
from scipy import linalg
from scipy.interpolate import RectBivariateSpline, interpn
import pdb

# Compute Focusing Delays
def calc_times(foci, elempos, dc = 0, speed_of_sound = 1540):
    ''' foc_times = calc_times(foci, elempos, dc = 0, speed_of_sound = 1540)

    CALC_TIMES - computes focusing times

    The function computes the (Tx or Rx) time of arrival for specified focal points
    given the array element positions.

    NOTE: Primarily intended when Tx and Rx apertures are the same (i.e. no full synthetic aperture)

    INPUTS:
    foci              - M x 3 matrix with position of focal points of interest [m]
    elempos           - N x 3 matrix with element positions [m]
    dc                - time offset [s]; scalar, N x 1 vector, or M x N array
    speed_of_sound    - speed of sounds [m/s]; default 1540 m/s

    OUTPUT:
    foc_times         - M x N matrix with times of flight for all foci and all array elements '''

    if type(dc).__module__ == 'builtins':
        dc = np.array([dc]);
    if not(np.isscalar(dc)) and sum(np.array(dc.shape)==1) <= 1:
        np.tile(dc, (foci.shape[0], 1));

    foci_tmp = np.tile(np.reshape(foci,(foci.shape[0],1,3)), (1,elempos.shape[0],1));
    elempos_tmp = np.tile(np.reshape(elempos,(1,elempos.shape[0],3)), (foci_tmp.shape[0],1,1));

    r = foci_tmp - elempos_tmp;

    distance = np.sqrt(np.sum(r**2, axis = 2));
    foc_times = distance/speed_of_sound + dc;

    return foc_times;


# Focus the RF Channel Data to Collect Receive Channel Data for Each Imaging Point
def tx_focus_fs(t, signal, foc_pts, rxAptPos, txAptPos = None, txApod = None, dc_rx = 0, dc_tx = 0, speed_of_sound = 1540):
    '''foc_data = tx_focus_fs(t, signal, foc_pts, rxAptPos, txAptPos = None, txApod = None, dc_rx = 0, dc_tx = 0, speed_of_sound = 1540)

    TX_FOCUS_FS - Focused the RF data at desired locations 
    (Full Delay-and-Sum on Transmit; Delayed but not Summed on Receive)

    The function interpolates the RF signals collected using the full synthetic sequence
    to focus the data at desired locations.

    INPUTS:
    t                  - T x 1 time vector for samples of the input signal
    signal             - T x N x M matrix containing input RF data to be interpolated
    foc_pts            - P x 3 matrix with position of focal points [m]
    rxAptPos           - N x 3 matrix with positions of the Rx apertures (elements) [m]
    txAptPos           - M x 3 matrix with positions of the Tx apertures (elements) [m]
                       - txAptPos = rxAptPos by default
    txApod             - P x M matrix of transmit apodizations for each Rx element and focal point
    dc_rx, dc_tx       - time offsets [s] for Tx and Rx; scalars, N (M) x 1 vectors, or P x N (M) matrix
    speed_of_sound     - speed of sounds [m/s]; default 1540 m/s

    OUTPUT:
    foc_data - vector with dimension P for beamformed image '''

    # Set variables to defaults if not set
    if txApod is None: txApod = np.ones((foc_pts.shape[0], rxAptPos.shape[0]));
    if txAptPos is None: txAptPos = rxAptPos;

    # time from the focus to receive  apertures (array elements)
    rx_times = calc_times(foc_pts, rxAptPos, dc = dc_rx, speed_of_sound = speed_of_sound);

    # time from the transmit apertures (array elements) to focus
    tx_times = calc_times(foc_pts, txAptPos, dc = dc_tx, speed_of_sound = speed_of_sound);

    # focused but not summed rf data
    foc_data = np.zeros((foc_pts.shape[0],rxAptPos.shape[0])).astype('complex64');
    for i in np.arange(rx_times.shape[1]):
        for j in np.arange(tx_times.shape[1]):
            foc_data[:,i] = foc_data[:,i] + txApod[:,j] * \
                np.interp(rx_times[:,i]+tx_times[:,j], t, signal[:,i,j], left=0, right=0);
    return foc_data;


# Generate Covariance Matrices for MIST
def MISTcov(nRx):
    '''mainlobe_cov, sidelobe_cov, noise_cov = MISTcov(nRx)
    MISTcov - Generate Covariance Matrices for Mainlobe, Sidelobes, and Noise
    Inputs:
        nRx = number of recieve channels
    Outputs:
        mainlobe_cov = covariance matrix for mainlobe
        sidelobe_cov = covariance matrix for sidelobes
        noise_cov = covariance matrix for noise '''
    
    # Sampling Rate (Samples Per Element)
    Fs = 10000; # Make this a Large Number for Pseudoanalytic Solution

    # Spatial and Wavenumber Axes
    Nx = Fs*nRx+1; # Total Number of X Samples
    kx = Fs*np.arange(-(Nx-1)/2,1+(Nx-1)/2)/Nx; # Arbitrary Lateral Coordinate for Far-Field Pattern

    # Create Wavenumber-Domain Signals and Their Inverse Fourier Transforms
    psf_vcz_orig = np.maximum(0,1-np.abs(kx)); # Linearly Decreasing VCZ curve
    box_vcz_orig = 2*np.sinc(2*kx); # Fourier Transform of Rect Function over the Mainlobe
    mainlobe_vcz_orig = np.real(np.fft.ifft(np.fft.fft(psf_vcz_orig, n=2*(Nx-1)+1)* 
        np.fft.fft(box_vcz_orig, n=2*(Nx-1)+1)))/nRx; # Convolve to Get Mainlobe Only
    mainlobe_vcz_orig = mainlobe_vcz_orig[int(Nx-1-(Nx-1)/2):int(Nx+(Nx-1)/2)];

    # Generate VCZ Curves for Mainlobe and Sidelobes 
    xcorr_aperture = np.logical_and(kx>=0, kx<=1);
    lag = kx[xcorr_aperture];
    full_vcz = np.real(psf_vcz_orig[xcorr_aperture]);
    mnlobe_vcz = np.real(mainlobe_vcz_orig[xcorr_aperture]);
    sdlobe_vcz = full_vcz-mnlobe_vcz;

    # Form Covariance Matrices for Mainlobe, Sidelobes, and Noise
    sidelobe_vcz = np.interp(np.arange(nRx)/nRx, lag, sdlobe_vcz);
    sidelobe_cov = linalg.toeplitz(sidelobe_vcz, sidelobe_vcz);
    mainlobe_vcz = np.interp(np.arange(nRx)/nRx, lag, mnlobe_vcz);
    mainlobe_cov = linalg.toeplitz(mainlobe_vcz, mainlobe_vcz);
    noise_cov = np.eye(nRx); # Covariance Matrix for Noise
    return mainlobe_cov, sidelobe_cov, noise_cov;


# Define Loadmat Function for HDF5 Format ('-v7.3' in MATLAB)
import h5py
def loadmat_hdf5(filename):
    file = h5py.File(filename,'r')
    out_dict = {}
    for key in file.keys():
        out_dict[key] = np.ndarray.transpose(np.array(file[key]));
    file.close()
    return out_dict;


# Python-Equivalent Command for IMAGESC in MATLAB
import matplotlib.pyplot as plt
def imagesc(x, y, img, rng, cmap='gray', numticks=(3, 3), aspect='equal'):
    exts = (np.min(x)-np.mean(np.diff(x)), np.max(x)+np.mean(np.diff(x)), \
        np.min(y)-np.mean(np.diff(y)), np.max(y)+np.mean(np.diff(y)));
    plt.imshow(np.flipud(img), cmap=cmap, extent=exts, vmin=rng[0], vmax=rng[1], aspect=aspect);
    plt.xticks(np.linspace(np.min(x), np.max(x), numticks[0]));
    plt.yticks(np.linspace(np.min(y), np.max(y), numticks[1]));
    plt.gca().invert_yaxis();
