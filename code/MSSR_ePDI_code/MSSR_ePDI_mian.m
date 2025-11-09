clc;
clear;
close all;


%--------------------------------------------------------------------------
% <<Enhanced Ultrasound Power Doppler Imaging via Mean-Shift Super-Resolution (MSSR)>>
% mian code
% Author           : Liu Y
% Date             : [2025-10-01]
% This code ia modified for ultrasound Power Doppler imaging from publicly accessible MSSR code at   https://github.com/adanog/MSSR/blob/main/MSSR_2.0.0_matlab.rar
% Reference        : 
% 1) Torres-García, E., Pinto-Cámara, R., Linares, A. et al. Extending resolution within a single imaging frame. Nat Commun 13, 7452 (2022). https://doi.org/10.1038/s41467-022-34693-9
%
% Description:
% This script demonstrates the application of a Mean-Shift Super-Resolution (MSSR)
% algorithm to enhance ultrasound Power Doppler Imaging (PDI) by processing
% clutter-filtered blood flow data obtained from SVD-based clutter suppression.
%
% The MSSR method enhances spatial resolution by exploiting sub-pixel
% information and denoising, improving visualization of microvascular flow.
%
% Input Data:
% - IQ: Complex quadrature-demodulated, beamformed ultrasound data
% - PRSSinfo: Structure containing parameters for clutter filtering and signal processing
%     * SVDrank    : 2-element vector specifying lower and upper singular value thresholds for SVD clutter removal
%     * HPfC       : High-pass filter cutoff frequency (unused here)
%     * NEQ        : Noise-equivalent quantification flag (set to 0)
%     * rFrame     : Effective frame rate after coherent compounding (Hz)
%
% Intermediate Variables:
% - sIQ          : SVD-processed blood signal (complex)
% - amp          : MSSR upsampling factor (typical values: 4-6)
% - psf_x, psf_y : Estimated axial and lateral Point Spread Function (PSF) widths (FWHM, pixels)
%                  For contrast-enhanced data, measured from isolated microbubbles;
%                  For non-contrast data, measured from sparse speckle regions.
% - order        : MSSR enhancement order (usually 2 to 10)
% - interp       : Interpolation method within MSSR ('fourier' or 'bicubic')
% - intNorm      : Boolean flag for intensity normalization (recommended true for PDI)
% - mesh         : Grid compensation parameter
% - excOL        : Outlier exclusion flag (currently unused)
% - OutLiersTh   : Threshold for outlier removal
% - is_downsamp  : Downsampling flag to reduce memory consumption (0 or 1)
%
% Outputs:
% - mag_iMSSR    : MSSR-enhanced magnitude data (3D matrix)
% - PDI_output   : Final power Doppler image after MSSR processing
%
% Usage Notes:
% - Load preprocessed sIQ data (SVD filtered blood flow signal)
% - Adjust parameters (e.g., psf_x, psf_y) according to dataset characteristics
% - Visualize results with appropriate dynamic range scaling (dB)
%
%--------------------------------------------------------------------------

% % Uncomment and adapt the following to load your data and parameters


%%%%%%%%%%%%%%%%%%%%%load IQ data and perform SVD %%%%%%%%%%%%%%%%%%%%%%%%%

% load kindey_IQ.mat
% 
% PRSSinfo.SVDrank=[30 250];
% PRSSinfo.HPfC=25;
% PRSSinfo.NEQ=0;
% PRSSinfo.rFrame=1000;
% [handles.sIQ, handles.sIQHP, handles.sIQHHP, handles.eqNoise]=IQ2sIQ(IQ(1:end,:,:),PRSSinfo);
% sIQ=handles.sIQ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%show IQ and sIQ data per frame%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:1:100
%     figure(3)
%     imagesc(db(abs(IQ(:,:,i))).^(1));  %xlabel,ylabel,
%     colormap("gray")
%     colorbar
%     pause(0.1)
% end
% 
% for i=1:1:100
%     figure(4)
%     imagesc((abs(sIQ(:,:,i))).^(1));  %xlabel,ylabel,
%     colormap("gray")
%     colorbar
%     pause(0.1)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% load sIQ data and show PD image; if you use IQ data and
%%%%%%%%%%%%%%%%% perform SVD, please comment those code %%%%%%%%%%%%%%%%

% load liver_sIQ.mat
load kidney_sIQ.mat

PDI1=mean(abs(sIQ(1:end,:,:)).^2,3);
figure;
imagesc(10*log10(PDI1./max(PDI1(:))));
colormap("gray")
clim([-30 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%--------------------------------------------------------------------------%
amp = 4;         % Upsampling factor for MSSR
psf_x = 3;         % Axial PSF width (FWHM in pixels) - example: 3 for contrast-enhanced kidney data; 4 for labe-free liver data
psf_y = 4;        % Lateral PSF width (FWHM in pixels) - example: 4 for contrast-enhanced kidney data; 6 for labe-free liver data
order = 6;            % MSSR enhancement order (2-10 typical)
mesh = 1;             % Grid compensation enabled
interp = 'bicubic';   % Interpolation method ('bicubic' or 'fourier')
intNorm = true;       % Enable intensity normalization
excOL = true;         % Enable outlier removal
OutLiersTh = 0.1;     % Outlier threshold parameter
is_downsamp = 1;      % Enable downsampling to reduce memory footprint
% -------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Separate magnitude and phase components of sIQ
mag_sIQ = abs(sIQ);
phase_sIQ = angle(sIQ);

% Apply MSSR super-resolution enhancement on magnitude data
[mag_iMSSR, PDI_output] = tMSSR(mag_sIQ, size(sIQ, 3), amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh, is_downsamp);

% Reconstruct complex enhanced signal by restoring phase
imgResult = mag_iMSSR .* exp(1i * phase_sIQ);

% Visualize MSSR-enhanced Power Doppler image (linear and dB scale)
figure;
imagesc(10*log10(PDI_output ./ max(PDI_output(:))));
colormap('gray');
colorbar;
clim([-50 0]);
title('MSSR-enhanced PDI (dB Scale)');

% For comparison: Upsampled conventional PDI
PDI1_upsampled = abs(imresize(PDI1, amp, 'bicubic'));
figure;
imagesc(10*log10(PDI1_upsampled ./ max(PDI1_upsampled(:))));
colormap('gray');
colorbar;
clim([-30 0]);
title('Conventional PDI Upsampled (dB Scale)');



















