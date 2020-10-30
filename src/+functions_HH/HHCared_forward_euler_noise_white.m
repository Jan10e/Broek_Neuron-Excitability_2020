function [v, n] = HHCared_forward_euler_noise_white(v0, gCa, C, v_sn, I_app, v_init, n_init, t, sigma)
% Simulation of reduced HH+Ca+ model

dt = t(2)-t(1);
eta = zeros(1,length(t));

%% Parameters

% Reversal potentials
ENa = 55;
EK = -77;
ECa = 85;
El = -54.4;

% Conductances
gl = 0.3; %0.3 or 0.4
gNa = 120;
gK = 36;

I_pump = -17;


%% Euler solve info
v = v_init;
n = n_init;


%% Generate noise vector
for i = 1:length(t)-1
    eta(i+1) = sigma*sqrt(dt)*randn;
end


%% Integrate using  Forward Euler
for ii = 1:length(eta)-1
    
    %Euler method to find the next voltage value
    v(ii +1) = v(ii) + dt .* (1/C) * (-gNa * functions.minf60_gen(v(ii))^3 ...
        * (0.89 - 1.1 * n(ii)) * (v(ii)-ENa) - gK * n(ii)^4 ...
        * (v(ii) - EK) - gl * (v(ii) - El) - gCa * n(ii)^3 ...
        * (v(ii) - ECa) + I_pump + I_app) + eta(ii);
    n(ii +1) = n(ii) + dt .* ((1/functions.taun60_ab(v_sn)) ...
        * (functions.ninf60_genV0(v(ii), v0) - n(ii)));
      
end

end

