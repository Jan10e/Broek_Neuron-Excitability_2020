function [v_sn, w_sn, I_app_sn] = FHNnc_sn_coordinates(v0, options)
%Calculation of coordinates of saddle-node (SN) bifurcation
%   The saddle node is the situation where the V-nullcline and
%   w-nullcline are tangential, irrespective of the epsilon value,
%   i.e. where det(J) = 0

w_inf1 = @(v) 2./(1 + exp(-5 * (v)));

% get v_sn coordinates when nullclines are tangential
[v_sn, ~, exitflag] = fsolve(@(V) (1-V^2).*(-1)-(-1).*((10*exp(5*v0 ...
    - 5*V))/(exp(5*v0 - 5*V) + 1)^2), -1, options);

if exitflag>0 && v_sn<0 %VSN<0 as there is another SN point further away but not relevant
    w_sn = w_inf1(v_sn-v0);
    I_app_sn = -(v_sn - (1/3) * v_sn^3 - w_sn);
else
    w_sn = NaN;
    I_app_sn = NaN;
end

end

