function [iMSSR,PDI_output] = tMSSR(varargin)
    % tMSSR  Process 3D data or file, supports different lateral and axial PSF parameters
    % Usage example:
    % [iMSSR, PDI_output] = tMSSR(inputData, dimz, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh)
    
    [inputData, dimz, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh, is_downsamp] = ParseInputs(varargin{:});

    % Initialize output variables
    if ischar(inputData) || isstring(inputData)
        % Input is a filename string, process frame by frame
        for j = 1:dimz
            Img = double(imread(inputData, j));
            iMSSR(:,:,j) = sfMSSR(Img, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh);
        end
        PDI_output = [];
    elseif isnumeric(inputData) && ndims(inputData) == 3
        % Input is a 3D numeric array, process directly
        [h,w,~] = size(inputData);
        PDI_large = zeros(h*amp, w*amp);

        for j = 1:dimz
             j
            Img = double(inputData(:,:,j));
            resAmp = sfMSSR(Img, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh);
            PDI_large = PDI_large + resAmp.^2;

            % Downscale back to original size
            if is_downsamp==1
                resSmall = imresize(resAmp, 1/amp, 'bicubic');
            else
                resSmall = resAmp;
            end

            iMSSR(:,:,j) = resSmall;

            clear Img resAmp resSmall
        end
        PDI_output = PDI_large ./ dimz;
    else
        error('tMSSR: The first input argument must be a filename string or a 3D numeric array!');
    end
end

function [inputData, dimz, amp, psf_x, psf_y, order, mesh, interp, intNorm, excOL, OutLiersTh, is_downsamp] = ParseInputs(varargin)
    narginchk(8,12); % Now requires at least 8 inputs, up to 12 inputs allowed

    inputData = varargin{1};
    dimz = varargin{2};
    amp = varargin{3};

    % Now psf_x and psf_y are the 4th and 5th input parameters respectively
    psf_x = varargin{4};
    psf_y = varargin{5};

    order = varargin{6};
    mesh = varargin{7};
    interp = varargin{8};
    interp = validatestring(interp, {'bicubic','fourier'}, mfilename, 'INTERP', 8);
    is_downsamp=varargin{12};

    % Set default values
    intNorm = false;
    excOL = false;
    OutLiersTh = 0.3;

    if length(varargin) >= 9
        intNorm = varargin{9};
    end
    if length(varargin) >= 10
        excOL = varargin{10};
    end
    if length(varargin) >= 11
        OutLiersTh = varargin{11};
    end
end





















% %tMSSR(imgName, dimz, amp, psf, order, mesh, interp, intNorm, excOL, OutLiersTh)
% function iMSSR = tMSSR(varargin)
%     disp(varargin);
%     [imgName, dimz, amp, psf, order, mesh, interp, intNorm, excOL, OutLiersTh] = ParseInputs(varargin{:});
%     for j = 1:dimz
% %         disp("Image " + j);
%         Img = double(imread(imgName,j));
%         [iMSSR(:,:,j)] = sfMSSR(Img, amp, psf, order, mesh, interp, intNorm, excOL, OutLiersTh);
%     end
% %     disp("--- Termine ---")
% end
% 
% function [imgName, dimz, amp, psf, order, mesh, interp, intNorm, excOL, OutLiersTh] = ParseInputs(varargin)
%     narginchk(7,10);
%     imgName = varargin{1};
%     dimz = varargin{2};
%     amp = varargin{3};
%     psf = varargin{4};
%     order = varargin{5};
%     mesh = varargin{6};
%     interp = varargin{7};
%     interp = validatestring(interp,{'bicubic','fourier'},mfilename,'INTERP',7);
%     if length(varargin) < 8
%         intNorm = false;
%         excOL = false;
%         OutLiersTh = 0.3;
%     elseif length(varargin) == 8
%         intNorm = varargin{8};
%         excOL = false;
%         OutLiersTh = 0.3;
%     elseif length(varargin) == 9
%         intNorm = varargin{8};
%         excOL = false;
%         OutLiersTh = 0.3;
%     elseif length(varargin) == 10
%         intNorm = varargin{8};
%         excOL = false;
%         OutLiersTh = varargin{10};
%     end
% end
