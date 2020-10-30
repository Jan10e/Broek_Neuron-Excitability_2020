function [v, n] = CSred_gen_forward_euler(gA, C, I_app, v_init, n_init, t)
%simulation of reduced Connor-Stevens model

%% Parameters

% Reversal potentials
ENa = 55;
EK = -75;
El = -17;

% Conductances
gl = 0.3;
gNa = 120;
gK = 20;


%% Euler solve info
dt = t(2)-t(1);

v = v_init;
n = n_init;


%% Solve ODE with Forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value
    v(i +1) = v(i) + dt .* (1/C).*(-gNa*functions.minf_gen(v(i))^3 ...
        * functions.hinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-ENa) ...
        - gK*n(i)^4*(v(i)-EK) - gA*functions.mAinf_gen(v(i))^3 ...
        * functions.hAinf_gen(functions.inv_ninf_gen(n(i)))*(v(i)-EK) ...
        - gl*(v(i)-El) + I_app);
    
    n(i +1) = n(i) + dt .* ((1/functions.tau_n(v(i))) ...
        * (functions.ninf_gen(v(i)) - n(i)));
    
end

end

