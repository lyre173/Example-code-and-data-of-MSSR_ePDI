%==========================================================================
% <<Enhanced Ultrasound Power Doppler Imaging via Mean-Shift Super-Resolution (MSSR)>>
% Ultrasound Flow Simulation in Single or Multiple Vessels with Anisotropic PSF
%
% Author           : Liu Y
% Date             : [2025-10-01]
%
%  - This script simulates blood flow inside single or multiple vessels
%    oriented horizontally or vertically, incorporating anisotropic PSF
%    (Point Spread Function) effects typical in ultrasound imaging.
%
% Description:
%  This MATLAB program generates simulated ultrasound sIQ data representing
%  moving scatterers within blood vessels. The vessels can be oriented
%  horizontally or vertically, with either single or double vessel modes.
%  The flow velocity profile is parabolic (centerline fast, edges slow),
%  and scatterers move accordingly frame-by-frame with wrap-around.
%  The imaging system's PSF is modeled as an anisotropic 2D Gaussian kernel
%  with different axial and lateral widths.
%
% Usage:
%  - Set 'mode_flag' to select vessel configuration:
%      1 - Single horizontal vessel
%      2 - Single vertical vessel
%      3 - Dual horizontal vessels
%      4 - Dual vertical vessels
%  - Configure parameters such as image size, PSF widths, vessel width,
%    vessel spacing, scatterers per vessel, maximum velocity, and frame count.
%  - Run the simulation to generate a 3D IQ data matrix representing moving
%    blood flow.
%  - Visualizes the flow simulation in real-time every 5 frames.
%
% Outputs:
%  - sIQ: Simulated IQ data cube (Nz x Nx x num_frames)
%
% Notes:
%  - Scatterer positions are updated according to velocity profiles and wrapped
%    at image boundaries to simulate continuous flow.
%  - PSF convolution models ultrasound imaging blur.
%  - Visualization uses jet colormap for intensity display.
%
% Requirements:
%  - MATLAB environment with basic image processing functions.
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% ========== Mode Selection ==========
% 1 - Single horizontal vessel
% 2 - Single vertical vessel
% 3 - Dual horizontal vessels
% 4 - Dual vertical vessels
mode_flag = 4;

%% Parameter Settings
Nx = 256;                  % Number of lateral pixels
Nz = 256;                  % Number of axial pixels
num_frames = 200;          % Total number of frames
sigma_x = 3;               % PSF lateral standard deviation (pixels)
sigma_z = 2;               % PSF axial standard deviation (pixels)
pipe_width = 3;            % Vessel width (pixels)
pipe_spacing = 6;          % Spacing between vessels (pixels)
num_points_per_pipe = 300; % Number of scatterers per vessel
v_max = 4;                 % Maximum velocity (pixels per frame)

%% Generate PSF (Anisotropic Gaussian)
[x_psf, z_psf] = meshgrid(-round(4*sigma_x):round(4*sigma_x), ...
                          -round(4*sigma_z):round(4*sigma_z));
PSF = exp(-(x_psf.^2/(2*sigma_x^2) + z_psf.^2/(2*sigma_z^2)));
PSF = PSF / sum(PSF(:)); % Normalize PSF energy

%% ========== Vessel Parameter Initialization ==========
switch mode_flag
    case 1 % Single horizontal vessel
        num_pipes = 1;
        pipe_centers = Nz/2;
        pipe_orient = 'horizontal';
    case 2 % Single vertical vessel
        num_pipes = 1;
        pipe_centers = Nx/2;
        pipe_orient = 'vertical';
    case 3 % Dual horizontal vessels
        num_pipes = 2;
        pipe_centers = [Nz/2 - pipe_spacing/2 - pipe_width/2, ...
                        Nz/2 + pipe_spacing/2 - pipe_width/2];
        pipe_orient = 'horizontal';
    case 4 % Dual vertical vessels
        num_pipes = 2;
        pipe_centers = [Nx/2 - pipe_spacing/2 - pipe_width/2, ...
                        Nx/2 + pipe_spacing/2 - pipe_width/2];
        pipe_orient = 'vertical';
    otherwise
        error('Invalid mode_flag, must be 1 to 4.');
end

%% ========== Initialize Scatterers Positions ==========
scatterers = struct();
for pipe = 1:num_pipes
    if strcmp(pipe_orient,'vertical')
        % Vertical vessel: x coordinate confined inside vessel width randomly,
        % z coordinate randomly distributed over entire depth
        x_initial = pipe_centers(pipe) - pipe_width/2 + pipe_width*rand(num_points_per_pipe,1);
        z_initial = Nz * rand(num_points_per_pipe,1);
    else
        % Horizontal vessel: z coordinate confined inside vessel width randomly,
        % x coordinate randomly distributed over entire lateral axis
        x_initial = Nx * rand(num_points_per_pipe,1);
        z_initial = pipe_centers(pipe) - pipe_width/2 + pipe_width*rand(num_points_per_pipe,1);
    end
    scatterers(pipe).x = x_initial;
    scatterers(pipe).z = z_initial;
end

%% ========== Visualization Initialization ==========
figure;
h_img = imagesc(zeros(Nz, Nx));
colormap(jet); axis image;
if strcmp(pipe_orient,'vertical')
    title(sprintf('%d Vertical Vessel(s) Flow Simulation', num_pipes));
    hold on;
    % Draw vessel boundaries as rectangles
    for i=1:num_pipes
        rectangle('Position',[pipe_centers(i)-pipe_width/2, 1, pipe_width, Nz], ...
                  'EdgeColor','w', 'LineStyle','--');
    end
else
    title(sprintf('%d Horizontal Vessel(s) Flow Simulation', num_pipes));
    hold on;
    for i=1:num_pipes
        rectangle('Position',[1, pipe_centers(i)-pipe_width/2, Nx, pipe_width], ...
                  'EdgeColor','w', 'LineStyle','--');
    end
end

%% ========== Main Loop ==========
IQ_frame = zeros(Nz,Nx,num_frames); % Pre-allocate IQ data cube

for frame = 1:num_frames
    current_scatter = zeros(Nz, Nx);

    for pipe = 1:num_pipes
        if strcmp(pipe_orient,'vertical')
            % Parabolic flow profile in vertical vessel:
            % faster at center, slower near edges
            x_offset = scatterers(pipe).x - pipe_centers(pipe);
            x_norm = x_offset / (pipe_width/2);
            velocities = v_max * (1 - x_norm.^2);
            
            % Update scatterer axial positions (z-direction)
            scatterers(pipe).z = scatterers(pipe).z + velocities;
            
            % Wrap around when exceeding bottom boundary
            overflow = scatterers(pipe).z > Nz;
            scatterers(pipe).z(overflow) = scatterers(pipe).z(overflow) - Nz;
            % Reinitialize lateral positions for overflowed scatterers inside vessel
            scatterers(pipe).x(overflow) = pipe_centers(pipe) - pipe_width/2 + ...
                                           pipe_width*rand(sum(overflow),1);
            % Convert to discrete pixel indices
            x_pos = round(scatterers(pipe).x);
            z_pos = mod(round(scatterers(pipe).z)-1, Nz)+1;
            valid = (x_pos >= 1) & (x_pos <= Nx);
        else
            % Parabolic flow profile in horizontal vessel
            z_offset = scatterers(pipe).z - pipe_centers(pipe);
            z_norm = z_offset / (pipe_width/2);
            velocities = v_max * (1 - z_norm.^2);
            
            % Update scatterer lateral positions (x-direction)
            scatterers(pipe).x = scatterers(pipe).x + velocities;
            
            % Wrap around when exceeding right boundary
            overflow = scatterers(pipe).x > Nx;
            scatterers(pipe).x(overflow) = scatterers(pipe).x(overflow) - Nx;
            % Reinitialize axial positions for overflowed scatterers inside vessel
            scatterers(pipe).z(overflow) = pipe_centers(pipe) - pipe_width/2 + ...
                                           pipe_width*rand(sum(overflow),1);
                                           
            % Convert to discrete pixel indices
            x_pos = mod(round(scatterers(pipe).x)-1, Nx)+1;
            z_pos = round(scatterers(pipe).z);
            valid = (z_pos >= 1) & (z_pos <= Nz);
        end
        % Keep only valid indices
        x_pos = x_pos(valid); 
        z_pos = z_pos(valid);
        
        % Convert subscript to linear indices
        lin_idx = sub2ind([Nz,Nx], z_pos, x_pos);
        
        % Accumulate scatterer counts in current frame
        current_scatter(lin_idx) = current_scatter(lin_idx) + 1;
    end
    
    % Convolve with PSF to simulate ultrasound imaging blur
    IQ_frame(:,:,frame) = conv2(current_scatter, PSF, 'same');
    
    % Real-time display every 5 frames
    if mod(frame,5) == 1
        set(h_img, 'CData', IQ_frame(:,:,frame));
        drawnow;
    end
end

% Assign output variable
sIQ = IQ_frame;

disp('Simulation completed!');














