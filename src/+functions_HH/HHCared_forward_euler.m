function [v, n] = HHCared_forward_euler(v0, gCa, C, v_sn, I_app, v_init, n_init, t)
% Simulation of reduced HH+Ca+ model

%% Constants

%Reversal potentials
ENa = 55;
EK = -77;
ECa = 85;
El = -54.4;

% conductances
gl = 0.3; %0.3 or 0.4
gNa = 120;
gK = 36;

Ipump = -17;

%% Euler solve info
dt = t(2)-t(1);

% initial values
% V = -70; %-40 or -70;
% % n = ninf60_genV0(V, V0);
% n = 0.2;

v = v_init;
n = n_init;

%% Integrate using  Forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value
    v(i +1) = v(i) + dt .* (1./C) * (-gNa * functions.minf60_gen(v(i))^3 ...
        * (0.89 - 1.1 * n(i)) * (v(i)-ENa) - gK * n(i)^4 ...
        * (v(i) - EK) - gl * (v(i) - El) - gCa * n(i)^3 ...
        * (v(i) - ECa) + Ipump + I_app);
    n(i +1) = n(i) + dt .* ((1/functions.taun60_ab(v_sn)) ...
        * (functions.ninf60_genV0(v(i), v0) - n(i)));
    
end

end

