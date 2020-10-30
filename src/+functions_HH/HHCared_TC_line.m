function gCa_TC = HHCared_TC_line(v0_range, options)
% Calculates transcritical bifurcation coordinates in 
% reduced HH+Ca

%% Constants

% Reversal potentials
ENa = 55;
EK = -77;
ECa = 85;
El = -54.4;

% Conductances
gl = 0.3;
gNa = 120;
gK = 36;

I_pump = -17;


%% Function

% Partials with V0 as dynamical variable and without timescale
dVdVC_V0 = @(V, n, gCa) (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - gCa*n^3 - 36*n^4 + (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) - 3/10;
dVdnC_V0 = @(V, n, gCa) (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - 144*n^3*(V + 77) - 3*gCa*n^2*(V - 85);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

% Partials of steady state values
dm3dV = @(V) (844424930131968*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4);
dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);

% gCa function
gCafun = @(V, n, V0) (1/n^3) * (-gNa*(0.89-1.1*n) ...
    * (dm3dV(V) * (V-ENa) + functions_HH.minf60_gen(V)^3) ...
    - gK * n^4 - gl);

%% loop gCa-TC
for i = 1:length(v0_range)
    v0 = v0_range(i);
    
    %get coordinates of TC using the expression for gAfun
    v_TC = fsolve(@(V) dVdnC_V0(V, functions_HH.ninf60_genV0(V, v0), ...
        gCafun(V, functions_HH.ninf60_genV0(V, v0), v0)) ...
        * dninfdV_V0(V, v0), -53, options);
    n_TC = functions_HH.ninf60_genV0(v_TC, v0);

    %get gA and Iapp at TC
    gCa_TC(i) = fsolve(@(gCa) dVdVC_V0(v_TC, n_TC, gCa), 0, options);

end

end

