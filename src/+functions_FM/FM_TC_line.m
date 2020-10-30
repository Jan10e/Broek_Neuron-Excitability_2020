function w0_TC = FM_TC_line(v0_range, options)
%Calculates transcritical bifurcation coordinates in mFHN model
%   TC values as mentioned in Franci(Siam, 2012)

I_app_TC = 2/3;
v_TC = -1;

for i = 1:length(v0_range)
    v0 = v0_range(i);
    
    w0_TC(i) = fsolve(@(w0) v_TC - (v_TC^3/3) ...
        - (2./(1 + exp(-5.* (v_TC - v0))) + w0 )^2 ...
        + I_app_TC, 0, options);
end

end

