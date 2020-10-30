function [v, m, h, n, d] = HHCafull_forward_euler_60mV(C, I_app, v_init, t)
%Simulation full model of HH+Ca
%   Use the 60mV adjusted values

%% Parameters

% Reversal potentials
ENa = 55;       %(shifted -60mV)
EK = -77;       %(shifted -60mV)
El = -54.4;     %(shifted -60mV)
ECa = 85;       %(shifted -60mV)

% ENa = 120;    %original HH
% EK = -12;     %original HH
% El = 10.6;    %original HH
% ECa = 150;

% Conductances  mS/cm^2
gNa = 120;
gK = 36;
gl = 0.3;
gCa = 0.4;

Ipump = -17;    % see Drion(2012)

%% Euler solve info
dt = t(2)-t(1);

%initial values
v = v_init;
m = functions.minf60_gen(v);   % Sodium current activation variable
h = functions.hinf60_gen(v);   % Sodium current inactivation variable
n = functions.ninf60_gen(v);   % Potassium current activation variables
d = functions.dinf60_gen(v);   % Calcium current activation variables


%% Integrate with forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next gating variable value
    m(i+1) = m(i) + dt .* ((1/functions.taum60_ab(v(i)))*(functions.minf60_gen(v(i)) - m(i)));
    n(i+1) = n(i) + dt .*((1/functions.taun60_ab(v(i)))*(functions.ninf60_gen(v(i)) - n(i)));
    h(i+1) = h(i) + dt .*((1/functions.tauh60_ab(v(i)))*(functions.hinf60_gen(v(i)) - h(i)));
    d(i+1) = d(i) + dt .*((1/functions.taud60(v(i)))*(functions.dinf60_gen(v(i)) - d(i)));
    
    %Euler method to find the next voltage value
    v(i+1) = v(i) + dt .*(1/C).*(-gNa*m(i)^3*h(i)*(v(i)-ENa) ...
        -gK*n(i)^4*(v(i)-EK) -gl*(v(i)-El) -gCa*d(i)*(v(i)-ECa) + Ipump + I_app );
    
end

end

