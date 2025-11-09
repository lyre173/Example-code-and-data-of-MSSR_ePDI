%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean Shift Core (Elliptical Kernel, Supports Lateral and Axial PSF) %%%%%%%%%%%%%%%%%%%%%%%%%%%

function MS = ms_core(I, psf_x, psf_y, amp)
% ms_core_ultrasound - Mean shift core function with elliptical kernel for ultrasound imaging
%
% Inputs:
%   I      - Input 2D image matrix (ultrasound image)
%   psf_x  - Lateral point spread function width (pixels)
%   psf_y  - Axial point spread function width (pixels)
%   amp    - Amplification factor, adjusts neighborhood size
%
% Output:
%   MS     - Processed mean shift result image
%
% Description:
%   This function accounts for the anisotropic PSF in ultrasound imaging by
%   applying an elliptical spatial kernel combined with an intensity kernel
%   for weighted mean shift calculation. The lateral and axial neighborhood
%   weights are calculated independently. The combination of spatial and 
%   intensity weights enhances edges and target resolution.

% Calculate anisotropic neighborhood radii (minimum 1 pixel)
hs_x = max(round(0.5 * psf_x * amp), 1);
hs_y = max(round(0.5 * psf_y * amp), 1);

% Use maximum radius for symmetric padding to ensure boundary safety
hs_max = max(hs_x, hs_y);
xPad = padarray(I, [hs_max, hs_max], 'symmetric');

[height, width] = size(I);

% Define neighborhood index ranges (rows = axial, columns = lateral)
int_y = -hs_y:hs_y;  % Axial range
int_x = -hs_x:hs_x;  % Lateral range

% Initialize matrix for maximum intensity difference calculation (for range normalization)
M = zeros(size(I));

% ------- First pass: compute max absolute intensity difference in neighborhood -------
% Purpose: normalize intensity difference based on local max difference to avoid scale inconsistency in weights
for dy = int_y
    for dx = int_x
        % Elliptical neighborhood condition (normalized squared sum ≤ 1)
        if (dy ~= 0 || dx ~= 0) && ((dy/hs_y)^2 + (dx/hs_x)^2 <= 1)
            % Extract neighborhood patch
            patch = xPad(hs_max + (1:height) + dy, hs_max + (1:width) + dx);
            % Update max difference matrix
            M = max(M, abs(I - patch));
        end
    end
end

% Avoid divide-by-zero by setting minimum difference threshold
M(M < 1e-6) = 1e-6;

% Initialize accumulators for weights and weighted intensity sum
weightAccum = zeros(size(I));
yAccum = zeros(size(I));

% ------- Second pass: compute weighted mean shift -------
% Combine spatial kernel (elliptical Gaussian) and range kernel (intensity difference based)
for dy = int_y
    for dx = int_x
        if (dy ~= 0 || dx ~= 0) && ((dy/hs_y)^2 + (dx/hs_x)^2 <= 1)
            % Compute spatial kernel weight (elliptical Gaussian)
            spatialKernel = exp(-0.5 * ((dy/hs_y)^2 + (dx/hs_x)^2));

            % Extract corresponding neighborhood patch
            patch = xPad(hs_max + (1:height) + dy, hs_max + (1:width) + dx);

            % Compute squared normalized intensity difference
            diffNormSq = ((I - patch) ./ M).^2;

            % Range kernel weight based on intensity difference (Gaussian)
            intensityKernel = exp(-0.5 * diffNormSq);

            % Combine weights
            weight = spatialKernel .* intensityKernel;

            % Accumulate weighted sums
            weightAccum = weightAccum + weight;
            yAccum = yAccum + patch .* weight;
        end
    end
end

% Compute mean shift vector (original image minus weighted mean)
MS = I - (yAccum ./ weightAccum);

% Numerical corrections: set negative values to zero to suppress noise
MS(MS < 0) = 0;

% Replace possible NaNs (from division by zero) with zero
MS(isnan(MS)) = 0;

end



