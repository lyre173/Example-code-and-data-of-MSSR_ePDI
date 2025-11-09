%--------------------------------------------------------------------------
% <<Enhanced Ultrasound Power Doppler Imaging via Mean-Shift Super-Resolution (MSSR)>>
% Functional Ultrasound Power Doppler Time-Series Analysis using Sliding Window
% mian code
% Author           : Liu Y
% Date             : [2025-10-01]
%
% Description:
% This script performs time-resolved analysis of functional ultrasound power 
% Doppler imaging (PDI) using a sliding window approach. It processes 
% beamformed IQ data files sequentially, applies SVD-based clutter suppression, 
% and optionally enhances the power Doppler signal with Mean-Shift Super-Resolution 
% (MSSR) processing (MSSR_ePDI).
%
% Main Features:
% - Batch loading of sequential IQ datasets with configurable file naming
% - Spatial cropping and temporal frame selection
% - Singular Value Decomposition (SVD) clutter filtering with user-defined thresholds
% - Interactive ROI (Region of Interest) selection on the first file's MSSR-ePDI
% - Sliding window computations of mean PDI intensity inside ROI over time
% - Comparison of conventional and MSSR-ePDI intensity time courses
% - Visualization of results and saving outputs as MAT and CSV files
% - Statistical correlation analysis between conventional and MSSR time series
%
% Input Parameters:
% - dataDir    : Directory path containing IQ .mat files
% - filePrefix : Prefix string for filenames
% - fileDigits : Number of digits in file numbering (e.g., 3 for "001")
% - fileStart  : Starting file index
% - fileEnd    : Ending file index
% - fileExt    : File extension (typically '.mat')
%
% - PRSSinfo: Structure with clutter suppression parameters
%    * SVDrank : Lower and upper SVD rank thresholds for clutter removal
%    * HPfC    : High-pass filter cutoff frequency (unused here)
%    * NEQ     : Noise-equivalent quantification flag (set to 0)
%    * rFrame  : Effective frame rate after compounding (Hz)
%
% - ROI cropping parameters (rows and columns)
% - useFrames : Number of frames used per file for analysis
%
% - Sliding window parameters:
%    * winSize  : Window length in frames
%    * stepSize : Sliding step size in frames
%
% - MSSR super-resolution parameters:
%    * amp        : Upsampling factor
%    * psf_x      : PSF axial full-width at half maximum (FWHM) in pixels
%    * psf_y      : PSF lateral FWHM in pixels
%    * order      : MSSR enhancement order (typically 2-10)
%    * mesh       : Grid compensation flag
%    * interp     : Interpolation method ('bicubic' or 'fourier')
%    * intNorm    : Intensity normalization flag
%    * excOL      : Outlier exclusion flag
%    * OutLiersTh : Outlier threshold
%    * is_downsamp: Downsampling flag for memory optimization (0 or 1)
%
% Outputs:
% - Time-course of mean PDI intensity within user-defined ROI for both conventional 
%   and MSSR-enhanced images
% - Figures displaying intensity vs time curves, ROI selection on PDI image
% - MAT and CSV files containing results and metadata
% - Statistical Pearson correlation coefficient and 95% confidence interval
%
% Usage Notes:
% - Requires IQ data pre-stored in .mat files with variable 'IQ' [rows x cols x frames]
% - Interactive ROI selection appears on the first file's MSSR-processed average PDI
%
% See also: IQ2sIQ, tMSSR
%--------------------------------------------------------------------------





clear; clc; close all;

%================ User Configuration ====================
dataDir     = 'E:\PDI20250926\IQ10'; % Directory containing IQ files
filePrefix  = 'liver_';              % Filename prefix
fileDigits  = 3;                     % Number of digits in filename numbering (e.g., '001')
fileStart   = 1;                     % Starting file number
fileEnd     = 52;                    % Ending file number
fileExt     = '.mat';                % File extension

% SVD clutter suppression parameters
PRSSinfo.SVDrank = [30 350];  
PRSSinfo.HPfC    = 25;    % Not used in this script
PRSSinfo.NEQ     = 0;
PRSSinfo.rFrame  = 1000;  % Frame rate after compounding (Hz)

% Spatial cropping parameters
rowCropStart = 15;           % Crop start row
rowCropEndPad= 70;           % Crop end from bottom (pad)
colCropStart = 25;           % Crop start column
colCropEnd   = 80;           % Crop end column (130-50)

useFrames    = 800;          % Number of frames used per file

% Sliding window parameters
winSize   = 200;             % Window length in frames
stepSize  = 20;              % Sliding step size in frames

% Upscaling factor for image resizing (for visualization)
re_size = 4;                 

% MSSR super-resolution parameters
amp = 4;                    % Upsampling factor
psf_x = 4;                  % Axial PSF width (FWHM, pixels)
psf_y = 3;                  % Lateral PSF width (FWHM, pixels)
order = 6;                  % MSSR enhancement order
mesh = 1;                   % Grid compensation enabled
interp = 'bicubic';         % Interpolation method ('bicubic' or 'fourier')
intNorm = true;             % Intensity normalization enabled
excOL = true;               % Outlier exclusion enabled
OutLiersTh = 0.3;           % Outlier exclusion threshold
is_downsamp = 0;            % Disable downsampling for MSSR

% Initialize results storage variables
intensityVals = [];
intensityVals_mssr = [];
timeVals = [];

% Initialize ROI mask variables
maskROI = [];
maskRefSize = [];

% Global frame counter for continuous timing across files
globalFrameCount = 0;

%================ Processing files sequentially =========================
for fi = fileStart:fileEnd
    % Construct the filename
    fname = fullfile(dataDir, [filePrefix, sprintf(['%0' num2str(fileDigits) 'd'], fi), fileExt]);
    fprintf('Processing file %d/%d : %s\n', fi - fileStart + 1, fileEnd - fileStart + 1, fname);
    if ~exist(fname, 'file')
        error('File does not exist: %s', fname);
    end
    
    % Load IQ data
    S = load(fname);
    if ~isfield(S, 'IQ')
        error('Variable ''IQ'' not found in file: %s', fname);
    end
    IQ = S.IQ;
    clear S;
    
    [nRow, nCol, nFr] = size(IQ);
    if nFr < useFrames
        error('File %s has insufficient frames: %d < %d', fname, nFr, useFrames);
    end
    
    % Crop spatial region of interest and select frames
    IQcrop = IQ(rowCropStart:(nRow - rowCropEndPad), colCropStart:colCropEnd, 1:useFrames);

    % Apply SVD-based clutter suppression to entire 800-frame block
    [handles.sIQ, handles.sIQHP, handles.sIQHHP, handles.eqNoise] = IQ2sIQ(IQcrop, PRSSinfo);
    sIQ = handles.sIQ;  % Clutter-filtered blood signal
    
    [sr, sc, sf] = size(sIQ);
    
    % Upscale conventional SVD-filtered images for visualization
    for t = 1:sf
        sIQ1(:,:,t) = imresize(sIQ(:,:,t), re_size, 'bicubic');
    end
    
    % Extract magnitude for MSSR processing
    mag_sIQ = abs(sIQ);

    % Perform MSSR super-resolution enhancement on magnitude data
    [mag_iMSSR, PDI_output] = tMSSR(mag_sIQ, sf, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh, is_downsamp);

    % On the first file, interactively select ROI on MSSR-enhanced mean PDI
    if fi == fileStart
        hFig = figure('Name','First File Mean MSSR-Enhanced PDI - Select ROI','NumberTitle','off');
        imagesc(20*log10(PDI_output ./ max(PDI_output(:))));
        colormap('gray');
        colorbar;
        clim([-60 0]);
        title('Draw an arbitrary ROI using mouse; Double-click or press Enter to finish');
        drawnow;
        
        % ROI tool selection depending on MATLAB version
        if exist('drawfreehand', 'file') == 2
            hROI = drawfreehand('LineWidth',1.5);
            maskROI = createMask(hROI);
        elseif exist('drawpolygon', 'file') == 2
            hROI = drawpolygon('LineWidth',1.5);
            maskROI = createMask(hROI);
        else
            maskROI = roipoly;
        end
        
        if isempty(maskROI)
            error('ROI mask not obtained. Please retry.');
        end
        maskROI = logical(maskROI);
        maskRefSize = size(maskROI);
        
        hold on; 
        visboundaries(maskROI, 'Color', 'w', 'LineWidth', 0.8); 
        hold off;
        fprintf('ROI saved and will be applied to all files.\n');
    end
    
    % Resize ROI mask if current image size differs
    if ~isequal([sr, sc], maskRefSize)
        maskResized = imresize(double(maskROI), [sr sc], 'nearest') > 0.5;
    else
        maskResized = maskROI;
    end
    
    % Sliding window computation of mean intensity inside ROI
    startFrames = 1:stepSize:(sf - winSize + 1);
    for startF = startFrames
        endF = startF + winSize - 1;
        
        % Conventional PDI block averaged intensity
        PDIblk = mean(abs(sIQ1(:,:,startF:endF)).^2, 3);
        ROIvals = PDIblk(maskResized);
        meanIntensity = mean(ROIvals(:));
        intensityVals(end+1, 1) = meanIntensity;

        % MSSR-enhanced PDI block averaged intensity
        PDIblk_mssr = mean(abs(mag_iMSSR(:,:,startF:endF)).^2, 3);
        ROIvals_mssr = PDIblk_mssr(maskResized);
        meanIntensity_mssr = mean(ROIvals_mssr(:));
        intensityVals_mssr(end+1, 1) = meanIntensity_mssr;

        % Compute time in seconds (assuming each useFrames = 1 second)
        timeVals(end+1, 1) = (globalFrameCount + startF + winSize/2) / useFrames; % Center time of window
    end
    
    globalFrameCount = globalFrameCount + sf;
    
    % Clear variables to save memory
    clear IQ IQcrop sIQ sIQ1 handles mag_iMSSR;
end

%================ Visualization of Time-Series ===========================
figure('Name','ROI Mean PDI Intensity vs Time (Sliding Window)','NumberTitle','off');
plot(timeVals, intensityVals, 'LineWidth', 1.4, 'MarkerSize', 4);
xlabel('Time (s)');
ylabel('Mean ROI PDI Intensity');
title('Functional Ultrasound - Power Doppler (Sliding Window)');
grid on;

figure('Name','ROI Mean MSSR-Enhanced PDI Intensity vs Time','NumberTitle','off');
plot(timeVals, intensityVals_mssr, 'LineWidth', 1.4, 'MarkerSize', 4);
xlabel('Time (s)');
ylabel('Mean ROI MSSR-Enhanced PDI Intensity');
title('Functional Ultrasound - MSSR-enhanced Power Doppler (Sliding Window)');
grid on;

% Save results to MAT and CSV files
outMat = fullfile(dataDir, 'PDI_timecourse_results_sliding.mat');
outCSV = fullfile(dataDir, 'PDI_timecourse_results_sliding.csv');
save(outMat, 'intensityVals', 'intensityVals_mssr', 'timeVals', 'maskROI', 'maskRefSize', 'PRSSinfo', 'winSize', 'stepSize');
T = table(timeVals, intensityVals, intensityVals_mssr);
writetable(T, outCSV);
fprintf('Results saved to:\n  %s\n  %s\n', outMat, outCSV);

% Save first file PDI + ROI figure
try
    saveas(hFig, fullfile(dataDir, 'PDI_firstfile_with_ROI.png'));
    fprintf('Saved first file PDI+ROI image.\n');
catch
    warning('Failed to save PDI+ROI image.');
end

%================ Relative Intensity Change Plot ==========================
figure;
plot(timeVals, 100 .* (intensityVals - mean(intensityVals(1:155))) ./ mean(intensityVals(1:155)), '-', 'LineWidth', 1.4, 'MarkerSize', 4);
hold on
plot(timeVals, 100 .* (intensityVals_mssr - mean(intensityVals_mssr(1:155))) ./ mean(intensityVals_mssr(1:155)), '-', 'LineWidth', 1.4, 'MarkerSize', 4);
xticks(0:5:52);
xlabel('Time [s]');
ylabel('Relative ROI PDI Change [%]');
title('Functional Ultrasound - MSSR-enhanced Power Doppler (Sliding Window)');
grid on;
xlim([0 52]);
legend('Conventional PDI', 'MSSR-enhanced PDI');

%================ Pearson Correlation and Confidence Interval =============
% Calculate relative percentage changes for correlation
PDI_data = 100 .* (intensityVals - mean(intensityVals(1:155))) ./ mean(intensityVals(1:155));
MSSR_PDI_data = 100 .* (intensityVals_mssr - mean(intensityVals_mssr(1:155))) ./ mean(intensityVals_mssr(1:155));

% Remove zero entries
idx = (PDI_data ~= 0) & (MSSR_PDI_data ~= 0);
PDI_data = PDI_data(idx);
MSSR_PDI_data = MSSR_PDI_data(idx);

% Compute Pearson correlation coefficient and p-value
[r, p] = corr(PDI_data, MSSR_PDI_data, 'type', 'Pearson');
n = length(PDI_data);

% Fisher z-transform for confidence interval
z = 0.5 * log((1 + r) / (1 - r));
SE = 1 / sqrt(n - 3);
zCI = [z - 1.96 * SE, z + 1.96 * SE];
rCI = (exp(2 * zCI) - 1) ./ (exp(2 * zCI) + 1);

fprintf('Pearson correlation coefficient r = %.3f, p = %.4f\n', r, p);
fprintf('95%% confidence interval: [%.3f, %.3f]\n', rCI(1), rCI(2));





