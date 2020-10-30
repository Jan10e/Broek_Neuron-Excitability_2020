function [v_sn, n_sn, I_app_sn, C_crit] = VSNCcrit(gCa, V0)
%Gets the VSN, IappSN and C_crit values for CSred

%% Parameters

%Reversal potentials
VNa = 55;
VK = -77;
VCa = 85;
Vl = -54.4;

% conductances
gl = 0.3;
gNa = 120;
gK = 36;

Ipump = -17;

%% differentials

%with V0 as dynamical variable and without timescale
dVdVC_V0 = @(V, n, gCa, V0) (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - gCa*n^3 - 36*n^4 + (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) - 3/10;
dVdnC_V0 = @(V, n, gCa, V0) (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - 144*n^3*(V + 77) - 3*gCa*n^2*(V - 85);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);


%% Get Iapp at SN bifn
options = optimoptions('fsolve','Display','none','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,'TolX',1e-14);

%Get SNIC bifurcation info by looking to range of IVs
IV_range = -80:-50;
for i = 1:length(IV_range)
    IV = IV_range(i);
    [v_sn_s(i), fval(i), exitflag(i)] = fsolve(@(V) ((dVdVC_V0(V, ninf60_genV0(V, V0), gCa, V0) * dndnT_V0) - (dVdnC_V0(V, ninf60_genV0(V, V0), gCa, V0) * dndVT_V0(V, V0))), IV, options);
end

if any(exitflag > 0)
    v_sn_n = v_sn_s(exitflag > 0); 
    v_sn_temp = uniquetol(v_sn_n, 1e-1);
    for ii = 1:length(v_sn_temp)
        n_sn_temp(ii) = ninf60_genV0(v_sn_temp(ii), V0);

        %% Check gs, which should be close to 0
        gs_temp(ii) = -(dVdnC_V0(v_sn_temp(ii), n_sn_temp(ii), gCa, V0) * dninfdV_V0(v_sn_temp(ii), V0));
        [gs_min, gs_idx] = min(abs(gs_temp));
        gs = gs_temp(gs_idx);

        %Collect the SN values
        v_sn = v_sn_temp(gs_idx);
        n_sn = n_sn_temp(gs_idx);
        I_app_sn = -( (-gNa*minf60_gen(v_sn).^3*(0.89-1.1*n_sn)*(v_sn-VNa) - gK*n_sn.^4*(v_sn-VK) - gl*(v_sn-Vl) - gCa* n_sn.^3* (v_sn-VCa) + Ipump) );
    end
else
    n_sn = nan;
    I_app_sn = nan;
    disp('no SN; 1 fixed point')
end

% Get C_crit
C_crit = fsolve(@(C) dVdVC_V0(v_sn, n_sn, gCa, V0) - (C/taun60_ab(v_sn)), 0, options);

end

