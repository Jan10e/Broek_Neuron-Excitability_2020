function [freq, v, n] = CSred_genV0_freq(gA, v0, v_sn, C, I_app, v_init, n_init, t)
%Calculate firing rate based on simulation and get the FI curve data

%% (1.) Integrate to get V-values
[v, n] = functions.CSred_genV0_forward_euler(v0(1), gA(1), C, v_sn(1), I_app(1), v_init, n_init, t);

%% (2.) Get peaks
[~, v_max_idx] = findpeaks(v, 'MinPeakProminence',5);
warning('off','signal:findpeaks:largeMinPeakHeight');
%      uncommend when checking peaks
%     figure;
%     findpeaks(V(1,:), 'MinPeakProminence',3, 'MinPeakHeight',10)

%% (3.) Get time of max peaks
t_max = t(v_max_idx);

%for no oscillations, its finding negative peaks
%     t_max = t_max(v_max > 0);

%% (4.) Get interspike interval
ISI = zeros(length(t_max),1);
if length(t_max) > 1
    for i = 1:length(t_max)-1
        ISI(i) = t_max(i+1) - t_max(i);
    end
    freq = 1000/mean(ISI); % for Hz
else
    freq = 0;
end
end

