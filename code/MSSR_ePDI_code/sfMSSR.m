function IMSSR = sfMSSR(varargin)
    % sfMSSR - Core function for ultrasound image enhancement based on mean shift
    % Input parameters:
    % x0          - Input 2D image matrix (double)
    % amp         - Upscaling factor
    % psf_x       - Lateral point spread function width
    % psf_y       - Axial point spread function width
    % order       - MSSR calculation order (number of iterations)
    % mesh        - Grid compensation switch (1 = enabled, 0 = disabled)
    % interp      - Interpolation method, 'bicubic' or 'fourier'
    % intNorm     - Whether to perform intensity normalization (optional, default false)
    % excOL       - Whether to exclude outliers (optional, default false)
    % OutLiersTh  - Outlier exclusion threshold (percentage, default 0.3)
    
    [x0, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh] = ParseInputs(varargin{:});
    
    type = "meanshift";
    x1 = double(x0);
    
    % Image upscaling
    if amp > 1
        if strcmpi(interp, "bicubic")
            AMP = imresize(x1, amp); % Bicubic interpolation upscaling
        else
            AMP = FourierMag(x1, amp); % Fourier domain interpolation upscaling
        end
        if mesh == 1
            AMP = compGrid(AMP, ceil(amp/2), 1);
        end
    else
        AMP = x1;
    end
    
    % Mean shift core analysis, passing lateral and axial PSF parameters
    MS = CoreAnalysisType(type, AMP, psf_x, psf_y, amp);
    
    % Normalization
    I3 = MS / max(MS(:));
    x3 = AMP / max(AMP(:));
    
    % MSSR iterative calculation
    for i = 1:order
        I4 = x3 - I3;          % Difference
        % The following commented code is alternative weighting factor calculation
        % diff=I4;
        % diff(diff<=1e-8)=0;
        % factor= (diff)./(x3+eps);
        % factor1=ones(size(factor));
        % factor1(factor>1)=0;  % Set threshold 0.8
        % factor=factor1;
        factor= (x3 - I3)./(x3+eps);
        factor(factor>0.8)=0;  % Set threshold 0.8
        factor( factor~=0)=1;
        
        I5 = max(I4(:)) - I4;  % Complement of difference
        I5 = I5 / max(I5(:));  % Normalization
        I6 = I5 .* I3.*factor; % Intensity weighting
        I7 = I6 / max(I6(:));  % Normalization
        x3 = I3;
        I3 = I7;
    end
    
    % Clean NaN values
    I3(isnan(I3)) = 0;
    
    IMSSR = I3;
    
    % Outlier exclusion
    if excOL
        th = (100 - OutLiersTh)/100;
        % The outlier exclusion specific implementation is commented out
        % Uncomment and adjust if needed
        % [f,x] = ecdf(AMP(:));
        % mnX = min(x(f>th));
        % AMP(AMP>mnX) = mnX;
        % [f,x] = ecdf(IMSSR(:));
        % mnX = min(x(f>th));
        % IMSSR(IMSSR>mnX) = mnX;
    end
    
    % Multiply back to upscaled image with intensity normalization
    if intNorm
        IMSSR = IMSSR .* AMP;
    end
end

% ------------- Helper Functions --------------

% Grid compensation function
function imgHVC = compGrid(img, desp, prp)
    [height, width, ~] = size(img);
    imgPad = padarray(img, [desp, desp], 'symmetric'); % Symmetric boundary padding
    imgVI = imgPad((1+desp):(height+desp), 1:width);
    imgVD = imgPad((1+desp):(height+desp), (2*desp + 1):(2*desp + width));
    imgHI = imgPad(1:height, (1+desp):(width+desp));
    imgHD = imgPad((2*desp + 1):(2*desp + height), (1+desp):(width+desp));
    imgHVC = (img + prp*(imgHD + imgHI + imgVD + imgVI)) / (1 + (4 * prp));
end

% Fourier domain interpolation upscaling function
function iFM = FourierMag(Img, mg)
    img = fft2(Img);
    sz = size(img);
    mdX = ceil(sz(1)/2);
    mdY = ceil(sz(2)/2);
    szF = round(sz * mg);
    fM = zeros(szF);
    lnX = length((mdX + 1):sz(1));
    lnY = length((mdY + 1):sz(2));
    img = szF(1)/sz(1) * szF(2)/sz(2) * img;

    fM(1:mdX, 1:mdY) = img(1:mdX, 1:mdY);
    fM(1:mdX, (szF(2)-lnY+1):szF(2)) = img(1:mdX, (1+mdY):sz(2));
    fM((szF(1)-lnX+1):szF(1), 1:mdY) = img((1+mdX):sz(1), 1:mdY);
    fM((szF(1)-lnX+1):szF(1), (szF(2)-lnY+1):szF(2)) = img((1+mdX):sz(1), (1+mdY):sz(2));
    iFM = ifft2(fM, 'symmetric');
end

% Input parsing function (updated to support psf_x and psf_y)
function [x0, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh] = ParseInputs(varargin)
    narginchk(7, 10);
    x0 = varargin{1};
    amp = varargin{2};
    psf_x = varargin{3};
    psf_y = varargin{4};
    order = varargin{5};
    mesh = varargin{6};
    interp = varargin{7};
    interp = validatestring(interp, {'bicubic', 'fourier'}, mfilename, 'INTERP', 7);

    % Set default values
    intNorm = false;
    excOL = false;
    OutLiersTh = 0.3;

    if length(varargin) >= 8
        intNorm = varargin{8};
    end
    if length(varargin) >= 9
        excOL = varargin{9};
    end
    if length(varargin) == 10
        OutLiersTh = varargin{10};
    end
end

% Example CoreAnalysisType function calling mean shift core (should be implemented and on path)
function MS = CoreAnalysisType(type, Img, psf_x, psf_y, amp)
    % Here simply demonstrates calling ms_core function
    if strcmp(type, "meanshift")
        MS = ms_core(Img, psf_x, psf_y, amp);
    else
        error('CoreAnalysisType: currently only supports type = "meanshift"');
    end
end
