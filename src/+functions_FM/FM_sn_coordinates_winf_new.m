function [v_sn, w_sn, I_app_sn] = FM_sn_coordinates_winf_new(v0, w0, options)
%Calculation of coordinates of saddle-node (SN) bifurcation
%   The saddle node is the situation where the V-nullcline and
%   w-nullcline are tangential, irrespective of the epsilon value

%% Functions
w_inf3 = @(v, v0, w0) 4./(1 + exp(-4 .* (v - (1/4) - v0))) + w0;
% w_inf3 = @(v, v0, w0) 4./(1 + exp(-4 .* (v - v0))) + w0;

%% get v_sn coordinates when nullclines are tangential, i.e. det(J)=0
[v_sn, ~, exitflag_v_sn]=fsolve(@(v) (1-v.^2).*(-1)-(-2.*w_inf3(v,v0,w0))...
    .*((16.*exp(4.*v0 - 4.*v + 1))./(exp(4.*v0 - 4.*v + 1) + 1).^2), -1.2, options);
% [v_sn, ~, exitflag_v_sn]=fsolve(@(v) (1-v.^2).*(-1)-(-2.*w_inf3(v,v0,w0))...
%     .*((16.*exp(4.*v0 - 4.*v))./(exp(4.*v0 - 4.*v) + 1).^2), -1.2, options);

%% wSN and Iapp values
if exitflag_v_sn>0
    w_sn = w_inf3(v_sn, v0, w0);
    I_app_sn = -(v_sn - (1/3) .* v_sn.^3 ...
        - (4./(1 + exp(-4 .*(v_sn - (1/4) - v0))) + w0).^2);
%     I_app_sn = -(v_sn - (1/3) .* v_sn.^3 ...
%         - (4./(1 + exp(-4 .*(v_sn - v0))) + w0).^2);
else
    w_sn = nan;
    I_app_sn = nan;
end

end