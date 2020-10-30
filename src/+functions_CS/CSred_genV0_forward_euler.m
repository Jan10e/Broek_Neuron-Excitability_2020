function [v, n] = CSred_genV0_forward_euler(v0, gA, C, v_sn, I_app, v_init, n_init, t)
%simulation of reduced Connor-Stevens model

%% Parameters

%Reversal potentials
ENa = 55;
EK = -75;
El = -17;

% conductances
gl = 0.3;
gNa = 120;
gKDR = 20;


%% Euler solve info
dt = t(2)-t(1);

%Initial values
v = v_init;
n = n_init;


%% Integrate using  Forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value
    v(i +1) = v(i) + dt .* (1./C).*(-gNa*functions.minf_gen(v(i))^3 ...
        * functions.hinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-ENa) ...
        - gKDR*n(i)^4*(v(i)-EK) - gA*functions.mAinf_gen(v(i))^3 ...
        * functions.hAinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-EK) ...
        - gl*(v(i)-El) + I_app);
    n(i +1) = n(i) + dt .* ((1./functions.tau_n(v_sn))...
        *(functions.ninf_genV0(v(i), v0) - n(i)));

end

end

