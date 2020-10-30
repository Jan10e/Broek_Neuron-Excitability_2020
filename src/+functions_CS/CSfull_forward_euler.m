function [v, m, h, n, a, b] = CSfull_forward_euler(C, I_app, v_init, t)
%Simulation of Connor-Stevens full model
%   gA is also a potassium channel, therefore EKDR is used as well

%% Parameters

% Reversal potentials
ENa = 55;
EK = -72;
EA = -75;
El = -17;

% Conductances (mS/mm^2)
gNa = 120; 
gK = 20;
gA = 47.7;
gl = 0.3;


%% Euler solve info
dt = t(2)-t(1);

%Initial values
v = v_init;                     % Initial Membrane voltage
m = functions.minf_ab(v);       % Sodium current activation variable
h = functions.hinf_ab(v);       % Sodium current inactivation variable
n = functions.ninf_ab(v);       % Potassium current activation variables
a = functions.a_inf(v);         % Potassium A current activation variables
b = functions.b_inf(v);         % Potassium A current activation variables


%% Integrate with forward Euler
for i = 1:length(t)-1
    
    %Euler method to find the next voltage value
    %eqn from Abbott pg 224    
    m(i+1) = m(i) + dt .* ((1/functions.tau_m(v(i)))*(functions.minf_ab(v(i)) - m(i)));
    h(i+1) = h(i) + dt .* ((1/functions.tau_h(v(i)))*(functions.hinf_ab(v(i)) - h(i)));
    n(i+1) = n(i) + dt .* ((1/functions.tau_n(v(i)))*(functions.ninf_ab(v(i)) - n(i)));
    a(i+1) = a(i) + dt .* ((1/functions.tau_a(v(i)))*(functions.a_inf(v(i)) - a(i)));
    b(i+1) = b(i) + dt .* ((1/functions.tau_b(v(i)))*(functions.b_inf(v(i)) - b(i)));
      
    %Euler method to find the next voltage value
    %eqn from Abbott pg 197 
    v(i+1) = v(i) + dt .*(1/C).*(-gNa*m(i)^3*h(i)*(v(i)-ENa) -gK*n(i)^4*(v(i)-EK) ...
        -gA*a(i)^3*b(i)*(v(i)-EA) -gl*(v(i)-El) + I_app );
   
end


end

