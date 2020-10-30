function [freq, v, w] = FM_freq(v0, w0, epsilon, I_app, v_init, w_init, t)
% Calculate firing rate based on simulation and get the FI curve data

%% (1.) Integrate to get v values
[v, w] = functions_FM.FM_forward_euler(v0, w0, epsilon, I_app, v_init, w_init, t);

%% (2.) Get peaks
[v_max, v_max_idx] = findpeaks(v, 'MinPeakProminence',0.05, 'MinPeakHeight', 0.1);
warning('off','signal:findpeaks:largeMinPeakHeight');
%   uncommend when checking peaks
%     figure;
%     findpeaks(V(1,:), 'MinPeakProminence',0.05, 'MinPeakHeight',-0.8)

%% (3.) Get time of max peaks
t_max = t(v_max_idx);

% Exclude negative peaks: for no oscillations, its finding negative peaks
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

