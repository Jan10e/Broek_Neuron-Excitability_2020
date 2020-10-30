function [v, m, h, n ,mA, hA, mCa] = CSCafull_forward_euler(C, I_app, v_init, t)
%Simulation of Connor-Stevens with Ca full model
%   gA is also a potassium channel, therefore EKDR is used as well

%% Parameters

% Reversal potentials
ENa = 55;
EKDR = -75;
ECa = 120;
El = -17;

% Conductances
gNa = 120;
gKDR = 20;
gCa = 0;
gA = 47.7;
gl = 0.3;

tau_mCa = 235/100;


%% ODE
dot_V = @(V, m, h, n, mA, hA, mCa) (1/C)*(-gNa*m^3*h*(V-ENa) ...
    -gKDR*n^4*(V-EKDR) - gA*mA^3*hA*(V-EKDR) - gCa*mCa^2*(V-ECa) ...
    - gl*(V-El) + I_app);


%% Euler solve info
dt = t(2)-t(1);

%Initial values
v = zeros(size(t)); % Voltage
v(1) = v_init;             % initial value
n=zeros(size(t));   % n: potassium activation gating variable
n(1) = 0.0;         
m=zeros(size(t));   % m: sodium activation gating variable
m(1) = 0.0;         
h=zeros(size(t));   % h: sodim inactivation gating variable
h(1) = 0.0;         
mA=zeros(size(t));   % A-current activation gating variable
mA(1) = 0.0;         
hA=zeros(size(t));   % A-current inactivation gating variable
hA(1) = 0.0;         
mCa = zeros(size(t));
mCa(1) = 0.0;       


%% Integrate with forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value   
    v(i + 1) = v(i) + dt .* dot_V(v(i), m(i), h(i), n(i), mA(i), hA(i), mCa(i));
    
    %Euler method to find the next gating variable value
    m(i + 1) = m(i) + dt .* functions.dm(v(i), m(i));
    h(i + 1) = h(i) + dt .* functions.dh(v(i), h(i));
    n(i + 1) = n(i) + dt .* functions.dn(v(i), n(i));
    mA(i + 1) = mA(i) + dt .* functions.dmA(v(i), mA(i));
    hA(i + 1) = hA(i) + dt .* functions.dhA(v(i), hA(i));
    mCa(i + 1) = mCa(i) + dt .* functions.dmCa(v(i), mCa(i));
    
end



end

