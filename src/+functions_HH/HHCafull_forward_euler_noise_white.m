function [v, m, h, n, d] = HHCafull_forward_euler_noise_white(C, I_app, v_init, t, sigma)
%Simulation full model of HH+Ca
%   Use the 60mV adjusted values

%% Constants

%reversal potentials
ENa = 55;       %(shifted -60mV)
EK = -77;       %(shifted -60mV)
El = -54.4;     %(shifted -60mV)
ECa = 85;       %(shifted -60mV)

% ENa = 120;    %original HH
% EK = -12;     %original HH
% El = 10.6;    %original HH
% ECa = 150;

%conductances  mS/cm^2
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

%% Generate noise vector
eta = zeros(1,length(t));

for i = 1:length(t)-1
    eta(i+1) = sigma*sqrt(dt)*randn(1);
end

%% Integrate with noise
for ii = 1:length(t)-1
    
    %Euler method to find the next m/n/h value
    m(ii+1) = m(ii) + dt .* ((1/functions.taum60_ab(v(ii)))*(functions.minf60_gen(v(ii)) - m(ii)));
    n(ii+1) = n(ii) + dt .*((1/functions.taun60_ab(v(ii)))*(functions.ninf60_gen(v(ii)) - n(ii)));
    h(ii+1) = h(ii) + dt .*((1/functions.tauh60_ab(v(ii)))*(functions.hinf60_gen(v(ii)) - h(ii)));
    d(ii+1) = d(ii) + dt .*((1/functions.taud60(v(ii)))*(functions.dinf60_gen(v(ii)) - d(ii)));
    
    %Euler method to find the next voltage value
    v(ii+1) = v(ii) + dt .*(1/C).*(-gNa*m(ii)^3*h(ii)*(v(ii)-ENa) ...
        -gK*n(ii)^4*(v(ii)-EK) -gl*(v(ii)-El) -gCa*d(ii)*(v(ii)-ECa) + Ipump + I_app ) + eta(ii);
    
end

end

