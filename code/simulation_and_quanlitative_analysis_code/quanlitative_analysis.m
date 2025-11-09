%--------------------------------------------------------------------------
% <<Enhanced Ultrasound Power Doppler Imaging via Mean-Shift Super-Resolution (MSSR)>>
% Quantitative Evaluation of Microvascular Power Doppler Images (dB Normalization)
% Author           : Liu Y
% Date             : [2025-10-01]
% Description:
% This main script quantitatively evaluates Power Doppler Images (PDI) and
% MSSR-enhanced Power Doppler Images (MSSR_ePDI) in the dB domain (normalized such that the maximum is 0 dB).
% All calculations are performed in the dB space for direct comparison.
%
% Steps:
%   1. Interpolate sIQ data, calculate conventional PDI and normalize to dB.
%   2. Set spatial resolutions (dx, dz) according to the dataset.
%   3. Visualize and compare both conventional and MSSR-enhanced PDI.
%   4. Replace NaN/Inf values with valid minimums.
%   5. Compute global and local contrast metrics.
%   6. Interactive ROI and blood/tissue region selection for CNR analysis.
%   7. Profile (line scan) analysis for spatial resolution (FWHM).
%
% INPUTS:
%   - sIQ: (M x N x T) complex data [should be loaded in workspace]
%   - PDI_output: MSSR-enhanced PDI (should be loaded/calculated)
%
% OUTPUTS:
%   - Quantitative metrics: global contrast, local contrast, CNR, FWHM
%   - Diagnostic figures
%
% USAGE NOTES:
%   - Requires Image Processing Toolbox for ROI and morphological operations.
%   - Adjust dx, dz according to actual experiment.
%   - Interactive ROI/line selection required.
%
%--------------------------------------------------------------------------


%%%% Interpolation %%%%%%%%%%%%%
% nt=size(sIQ,3);
% for t = 1:nt
%     sIQ1(:,:,t) = imresize(sIQ(:,:,t), 6, 'bicubic');
% end
% 
% 
% PDI1 = mean(abs(sIQ1).^2, 3);
% PDI_samp1 = 10.*log10(PDI1 ./ max(PDI1(:)));

% Interpolation and PDI calculation
nt = size(sIQ,3);             % Number of frames
scale = 6;                    % Interpolation scale factor
PDI_sum = 0;
for t = 1:nt
    % Interpolate each frame and accumulate squared magnitude
    PDI_sum = PDI_sum + abs(imresize(sIQ(:,:,t), scale, 'bicubic')).^2;
end
PDI_samp1 = 10 .* log10(PDI_sum/nt ./ max(PDI_sum(:)/nt)); % dB normalization



% ====== Set spatial resolution (modify as needed for your data) ======
% For simulated data (default example)


% %%%%%%%contrast-enhanced heart data %%%
% dx =0.113; % （mm）
% dz = 0.113; % （mm）

% %%%%%%contrast-enhanced mouse tumor%%%%%%
% dx =0.0183; % （mm）
% dz = 0.0221; % （mm）

% %%%%%%contrast-enhanced rat brain%%%%%%
% dx =0.0143; % （mm）
% dz = 0.0167; % （mm）

% %%%%%%label- free mouse liver and kidney%%%%%%
% dx =0.0856/6/1; % （mm）
% dz = 0.0856/6/1; % （mm）

% %%%%%%simulation data%%%%%%
% dx =1/6/1; % （mm）
% dz = 1/6/1; % （mm）


% %%%%%%%%%%%label-free rat brain%%%%%%%
dz=(P.xCoor(2)-P.xCoor(1))/6;
dx=(P.zCoor(2)-P.zCoor(1))/6;


[nx, nz] = size(PDI_samp1);
x = (0:nx-1) * dx;                               % Lateral coordinates (mm)
z = (-nz/2:nz/2) * dz;                           % Depth coordinates (mm)



PDI2 = 10.*log10(PDI_output./max(PDI_output(:)));
figure;
imagesc(z,x,imgaussfilt(PDI2, 1));
colormap("hot")
clim([-30 0])
%  axis equal
xlabel('Lateral [mm]')
ylabel('Depth [mm]')



% ====== Replace NaN and Inf with minimum valid value for each image ======
for matIdx = 1:2
    if matIdx==1
        M = PDI_samp1;
    else
        M = PDI2;
    end
    validVals = M(~isnan(M) & ~isinf(M));
    minVal = min(validVals);
    M(isnan(M) | isinf(M)) = minVal;
    if matIdx==1
        PDI_samp1 = M;
    else
        PDI2 = M;
    end
end



% Prepare data for next steps
data_list = {PDI_samp1, PDI2};
data_name = {'PDI', 'MSSR_ePDI'};
data_lin = {PDI1 ./ max(PDI1(:)), PDI_output ./ max(PDI_output(:))};

%% 1. Global Contrast and RMS Contrast Analysis (ROI selection)
disp('========= PART 1: Global and RMS Contrast Analysis ===========');
choice = questdlg('Would you like to manually select a ROI for global contrast analysis?', ...
    'Global Contrast ROI Selection', 'Yes','No (Whole Image)','No (Whole Image)');
manual_select = strcmp(choice, 'Yes');

if manual_select
    figure; imagesc(PDI_samp1, [-25 0]); axis image; axis equal; colormap('hot'); colorbar;
    title('Please draw a freehand ROI on conventional PDI, double-click or press Enter to finish');
    hroi_glob = drawfreehand('Color','g','LineWidth',2);
    wait(hroi_glob);
    globalMask = createMask(hroi_glob);   % Save global ROI mask
    close;
else
    globalMask = true(size(PDI_samp1));   % Use the whole image
end

for k = 1:2
    I = data_list{k};
    I_sel = I(globalMask);
    I_mean = mean(I_sel(:));
    std_all = std(I_sel(:));
    rms_c = sqrt(mean(I_sel(:).^2) - I_mean^2);
    fprintf('%s:\n  Global Contrast (Std, dB): %.4f, RMS Contrast (dB): %.4f\n', ...
        data_name{k}, std_all, rms_c);
end
disp('-----------------------------------------------------');
pause(1);

%% 2. ROI and CNR Analysis (linked to previous ROI)
disp('========= PART 2: ROI and CNR Analysis ===========');

if ~manual_select
    % Whole image previously used, must select a new ROI for CNR
    disp('NOTE: Please select a ROI for CNR analysis...');
    figure; imagesc(PDI_samp1, [-25 0]); axis image; colormap('hot'); colorbar;
    title('Draw a freehand ROI on conventional PDI, double-click to finish');
    hroi = drawfreehand('Color','g','LineWidth',2);
    wait(hroi);
    roiMask = createMask(hroi);
    close;
else
    % Ask if user wants to reuse or reselect ROI
    roi_choice = questdlg('Redraw ROI for CNR analysis?', ...
        'ROI Selection', 'Yes (Redraw)', 'No (Reuse global ROI)', 'No (Reuse global ROI)');
    roi_reselect = strcmp(roi_choice, 'Yes (Redraw)');
    if roi_reselect
        figure; imagesc(PDI_samp1, [-70 0]); axis image; colormap('hot'); colorbar;
        title('Draw a freehand ROI for CNR analysis, double-click to finish');
        hroi = drawfreehand('Color','g','LineWidth',2);
        wait(hroi);
        roiMask = createMask(hroi);
        close;
    else
        roiMask = globalMask;
    end
end

se_blood = strel('disk',3);         % Blood region morphology
se_blood_dilate = strel('disk',2);  % Blood region dilation
se_tissue = strel('disk',3);        % Tissue region erosion

for k = 1:2
    I_dB = data_list{k};                % Current image (dB)
    I_lin = data_lin{k};                % Current image (linear)
    ROI_vals_dB = I_dB(roiMask);
    ROI_vals_lin = I_lin(roiMask);

    % Threshold for blood region segmentation
    prompt = {['Please input blood region threshold (dB) for ' data_name{k} ':']};
    dlg_title = ['Blood Region Threshold - ' data_name{k}];
    defaultans = {num2str(mean(ROI_vals_dB))};
    answer = inputdlg(prompt, dlg_title, 1, defaultans);
    if isempty(answer), disp('User cancelled.'); return; end
    threshold = str2double(answer{1});

    % Blood region initial mask
    blood_mask = false(size(I_dB));
    blood_roi_idx = find(roiMask);
    is_blood = ROI_vals_dB > threshold;
    blood_mask(blood_roi_idx(is_blood)) = true;

    % Morphological processing
    blood_mask = imerode(blood_mask, se_blood);
    blood_mask = imdilate(blood_mask, se_blood);
    blood_mask = bwareaopen(blood_mask, 5);

    % Tissue region mask
    blood_dil = imdilate(blood_mask, se_blood_dilate);
    tissue_mask = ~blood_dil;
    tissue_mask = imerode(tissue_mask, se_tissue);
    tissue_mask = tissue_mask & roiMask;
    tissue_mask = bwareaopen(tissue_mask, 5);

    % Compute contrast and CNR in linear scale
    blood_vals = I_lin(blood_mask);
    tissue_vals = I_lin(tissue_mask);

    mean_blood = mean(blood_vals);
    mean_tissue = mean(tissue_vals);
    std_tissue = std(tissue_vals);

    contrast_db = 10*log10(mean_blood) - 10*log10(mean_tissue);
    CNR_lin = (mean_blood - mean_tissue) / std_tissue;
    CNR_db = 10*log10(CNR_lin);

    fprintf('%s: Local Contrast = %.2f dB, CNR = %.2f dB\n', ...
        data_name{k}, contrast_db, CNR_db);

    % Visualization
    figure;
    subplot(1,3,1); imagesc(I_dB, [-20 0]); axis image off; colormap('hot');
    hold on; contour(roiMask, [0.5 0.5], 'g', 'LineWidth', 2);
    title([data_name{k},' Original + ROI']); colorbar;
    subplot(1,3,2); imagesc(blood_mask); axis image off; title('Blood Mask');
    subplot(1,3,3); imagesc(tissue_mask); axis image off; title('Tissue Mask');
end

disp('-----------------------------------------------------');
pause(1);

%% 3. Profile Analysis (draw line, FWHM in dB space)
disp('========= PART 3: Profile (Line Scan) Analysis ===========');
figure; imagesc(z, x, PDI2); axis image; axis equal; colormap('gray'); colorbar;
clim([-30 0]);
xlabel('Lateral [mm]');
ylabel('Depth [mm]');
title('Draw a line across a vessel, double-click to finish');
hL = drawline('Color', 'w', 'LineWidth', 1.5);
wait(hL);
line_pos = hL.Position; % [z1 x1; z2 x2]

% Sampling: increase density for FWHM accuracy
N = round(pdist(line_pos*200));
z_line = linspace(line_pos(1,1), line_pos(2,1), N); % Lateral (mm)
x_line = linspace(line_pos(1,2), line_pos(2,2), N); % Axial (mm)

% Axes for interpolation
nz = size(PDI_samp1,2);
nx = size(PDI_samp1,1);
if length(z) ~= nz
    z_axis = linspace(min(z), max(z), nz);
else
    z_axis = z;
end
if length(x) ~= nx
    x_axis = linspace(min(x), max(x), nx);
else
    x_axis = x;
end

zq = interp1(z_axis, 1:nz, z_line, 'linear', 'extrap');
xq = interp1(x_axis, 1:nx, x_line, 'linear', 'extrap');
dist = sqrt((z_line - z_line(1)).^2 + (x_line - x_line(1)).^2);

% Profile extraction for both images
profile_curves = cell(1,2);
for k = 1:2
    I = data_list{k};
    profile_curves{k} = interp2(I, zq, xq, 'linear');
end

% Overlayed profile curves
figure;
plot(dist, profile_curves{1}, 'b-', 'LineWidth', 2); hold on;
plot(dist, profile_curves{2}, 'r--', 'LineWidth', 2);
xlabel('Distance [mm]');
ylabel('Intensity [dB]');
legend(data_name, 'Location', 'best');
title('Profile Comparison along Line');
grid on; set(gca,'FontSize',12);

% FWHM calculation and peak visualization for each profile
for k = 1:2
    prof = profile_curves{k};
    figure; plot(dist, prof, 'b-', 'LineWidth',2); grid on;
    xlabel('Distance [mm]'); ylabel('Amplitude (dB)');
    title(['Profile - ',data_name{k}]);
    hold on;
    [pks,locs_idx] = findpeaks(prof, 'MinPeakProminence', 3);
    locs_mm = dist(locs_idx);
    plot(locs_mm, pks, 'ro','MarkerFaceColor','r');
    fwhm_mm = [];
    for i = 1:length(pks)
        pk = pks(i);
        center_idx = locs_idx(i);
        half_level = pk - 3;
        % Left boundary
        left_idx = find(prof(1:center_idx)<half_level,1,'last');
        if isempty(left_idx), left_idx=1; end
        if left_idx<center_idx && prof(left_idx)<half_level && prof(left_idx+1)>half_level
            left_x = interp1(prof(left_idx:left_idx+1), dist(left_idx:left_idx+1), half_level);
        else
            left_x = dist(left_idx);
        end
        % Right boundary
        right_idx = find(prof(center_idx:end)<half_level,1,'first');
        if isempty(right_idx)
            right_idx = length(prof)-center_idx+1;
        end
        right_idx = right_idx+center_idx-1;
        if right_idx>center_idx && prof(right_idx-1)>half_level && prof(right_idx)<half_level
            right_x = interp1(prof(right_idx-1:right_idx), dist(right_idx-1:right_idx), half_level);
        else
            right_x = dist(right_idx);
        end
        width = right_x - left_x;
        fwhm_mm(i) = width;
        plot([left_x right_x],[half_level half_level],'k--','LineWidth',2);
        plot([locs_mm(i) locs_mm(i)],[min(prof),pk],'m:');
        text(locs_mm(i), pk, sprintf('%.2fmm', width), ...
            'Color','k','FontSize',8,'VerticalAlignment','bottom');
    end
    % Output results
    disp('--------------------------------------------');
    fprintf('【%s】Detected Vessel Peaks: %d\n', data_name{k}, length(pks));
    if isempty(pks)
        disp('No valid peaks detected!');
    else
        for i = 1:length(pks)
            fprintf('  Peak %d: Position=%.2fmm, FWHM(-3dB)=%.2f mm\n', i, locs_mm(i), fwhm_mm(i));
        end
        fprintf('All FWHMs (-3dB) (mm): ');
        fprintf('%.2f  ', fwhm_mm);
        fprintf('\n');
    end
    disp('--------------------------------------------');
end

















