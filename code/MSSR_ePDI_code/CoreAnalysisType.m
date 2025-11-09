% Selecction core process
function CoreImage = CoreAnalysisType(type, I, psf_x,psf_y, amp)
    if lower(type) == "meanshift"
        CoreImage = I;
        for i = 1:1
            CoreImage = ms_core(CoreImage, psf_x/i, psf_y, amp);
        end
    end
end
