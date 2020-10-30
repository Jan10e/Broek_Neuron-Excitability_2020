function [freq0, I0] = FM_freq_min(v0, w0, epsilon, I_app_range, v_init, w_init, t)
% Calculate firing rate based on simulation and get the FI curve data

for i = 1:length(I_app_range)
    I_app = I_app_range(i);
    
    % get frequency
    freq = functions.FM_freq(v0, w0, epsilon, I_app, v_init, w_init, t);
    
    % minimum frequency
    if freq > 0
        freq0 = freq;
        I0 = I_app;
        break
    else
        freq0 = 0;
        I0 = 0;
    end
    
end

