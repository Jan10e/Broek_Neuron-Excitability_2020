function [freq, v, n] = HHCared_freq(gCa, v0, v_sn, C, I_app, v_init, n_init, t)
%Calculate firing rate based on simulation and get the FI curve data

%% (1.) Integrate to get V-values
[v, n] = functions.HHCared_forward_euler(v0, gCa, C, ...
    v_sn, I_app, v_init, n_init, t);

%% (2.) Get peaks
[v_max, v_max_idx] = findpeaks(v, 'MinPeakProminence', 20, ...
    'MinPeakHeight',10, 'MinPeakDistance', 500);
%   uncommend when checking peaks
%     figure;
%     findpeaks(V(1,:), 'MinPeakProminence',20, 'MinPeakHeight',10, 'MinPeakDistance',500)

%% (3.) Get time of max peaks
t_max = t(v_max_idx);

%for no oscillations, its finding negative peaks
%     t_max_pos = t_max(v_max > 0);

%% (4.) Get interspike interval
ISI = zeros(length(t_max),1);
if length(t_max) > 1
    for i = 1:length(t_max)-1
        ISI(i) = t_max(i+1) - t_max(i);
    end
    freq = 1000/mean(ISI); %for Hz
else
    freq = 0;
end

end

