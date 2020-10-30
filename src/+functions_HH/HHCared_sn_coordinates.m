function [v_sn, n_sn, I_app_sn] = HHCared_sn_coordinates(gCa, v0, options)
%Calculation of coordinates of saddle-node (SN) bifurcation
%   The saddle node is the situation where the V-nullcline and
%   n-nullcline are tangential, irrespective of the epsilon value

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

%% Differentials

% With V0 as dynamical variable and without timescale
dVdVC_V0 = @(V, n, gCa, V0) (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - gCa*n^3 - 36*n^4 + (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) - 3/10;
dVdnC_V0 = @(V, n, gCa, V0) (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - 144*n^3*(V + 77) - 3*gCa*n^2*(V - 85);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);


%% Get v_sn coordiates

% Get SNIC bifurcation info by looking to range of IVs
[v_sn, ~, exitflag] = fsolve(@(V) ((dVdVC_V0(V, functions.ninf60_genV0(V, v0), ...
    gCa, v0) * dndnT_V0) - (dVdnC_V0(V, functions.ninf60_genV0(V, v0), gCa, v0) ...
    * dndVT_V0(V, v0))), -63, options);

% IV_range = -80:-50;
%
% for i = 1:length(IV_range)
%     IV = IV_range(i);
%     [v_sn_s(i), ~, exitflag(i)] = fsolve(@(V) ...
%         ((dVdVC_V0(V, functions.ninf60_genV0(V, V0),gCa, V0) ...
%         * dndnT_V0) - (dVdnC_V0(V, functions.ninf60_genV0(V, V0),gCa, V0) ...
%         * dndVT_V0(V, V0))), IV, options);
% end

%% n_SN and Iapp_SN
if any(exitflag > 0) && v_sn > -80
    n_sn = functions.ninf60_genV0(v_sn, v0);
    I_app_sn = -( (-gNa * functions.minf60_gen(v_sn).^3 ...
        * (0.89-1.1*n_sn) * (v_sn-ENa) - gK * n_sn.^4 ...
        * (v_sn-EK) - gl*(v_sn-El) - gCa * n_sn.^3 ...
        * (v_sn-ECa) + I_pump) );
    
    
    
    
    %     v_sn_n = v_sn_s(exitflag > 0);
    %     v_sn_temp = uniquetol(v_sn_n, 1e-1);
    %     for ii = 1:length(v_sn_temp)
    %         n_sn_temp(ii) = functions.ninf60_genV0(v_sn_temp(ii), V0);
    %
    %         %% Check gs, which should be close to 0
    %         gs_temp(ii) = -(dVdnC_V0(v_sn_temp(ii), ...
    %             n_sn_temp(ii), gCa, V0) ...
    %             * dninfdV_V0(v_sn_temp(ii), V0));
    %         [~, gs_idx] = min(abs(gs_temp));
    % %         gs = gs_temp(gs_idx);
    %
    %         %Collect the SN values
    %         v_sn = v_sn_temp(gs_idx);
    %         n_sn = n_sn_temp(gs_idx);
    %         I_app_sn = -( (-gNa*functions.minf60_gen(v_sn).^3 ...
    %             *(0.89-1.1*n_sn)*(v_sn-ENa) - gK*n_sn.^4 ...
    %             *(v_sn-EK) - gl*(v_sn-El) - gCa* n_sn.^3 ...
    %             *(v_sn-ECa) + Ipump) );
    
% end
else
    v_sn = nan;
    n_sn = nan;
    I_app_sn = nan;
end

end

