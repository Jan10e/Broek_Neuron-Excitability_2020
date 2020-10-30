function fixed_pts = HHCared_fixed_points_multi(steady_state, options)
%Multiple fixed points
%   To obtain all the fixed points, look at all places where the
%   V-nullcline crosses the n-nullcline

iter = -121:3:0; % make sure -50 is not present
fixed_pts_temp = nan(size(iter));

for i = 1:length(iter)
    x1 = iter(i);
    [fixed_pts_temp(i), ~, exitflag_fixed_pts(i)] = ...
        fsolve(steady_state, x1, options);
    fixed_pts_temp_converged = exitflag_fixed_pts > 0;
    fixed_pts_temp_converged = fixed_pts_temp(fixed_pts_temp_converged);
end
fixed_pts = uniquetol(fixed_pts_temp_converged, 1E-4);

end

