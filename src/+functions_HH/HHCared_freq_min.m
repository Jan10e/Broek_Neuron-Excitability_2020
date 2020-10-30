function [freq0, I0] = HHCared_freq_min(gCa, v0, v_sn, C_crit, I_app_range, v_init, n_init, t)
%Calculate firing rate based on simulation and get the FI curve data

for i = 1:length(I_app_range)
    I_app = I_app_range(i);
    
    % get frequency
    freq = functions.HHCared_freq(gCa, v0, v_sn, C_crit, I_app, v_init, n_init, t);
    
    %minimum frequency
    if freq > 0
        freq0 = freq;
        I0 = I_app;
        break
    else
        freq0 = 0;
        I0 = 0;
    end

end

