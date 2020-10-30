function [v, n] = CSred_genV0_forward_euler_noise_white(v0, gA, C, v_sn, I_app, v_init, n_init, t, sigma)
%simulation of reduced Connor-Stevens model

%% Parameters

%Reversal potentials
ENa = 55;
EK = -75;
El = -17;

% conductances
gl = 0.3;
gNa = 120;
gK = 20;


%% Euler solve info
dt = t(2)-t(1);

% initial values
v = v_init;
n = n_init;

%% Generate noise vector
eta = zeros(1,length(t));

for i = 1:length(t)-1
    eta(i+1) = sigma*sqrt(dt)*randn;
end

%% Integrate using  Forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value
    v(i +1) = v(i) + dt .* (1./C).*(-gNa*functions.minf_gen(v(i))^3 ...
        * functions.hinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-ENa) ...
        - gK*n(i)^4*(v(i)-EK) - gA*functions.mAinf_gen(v(i))^3 ...
        * functions.hAinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-EK) ...
        - gl*(v(i)-El) + I_app) + eta(i);
    n(i +1) = n(i) + dt .* ((1./functions.tau_n(v_sn)) ...
        * (functions.ninf_genV0(v(i), v0) - n(i)));
    
end
end

