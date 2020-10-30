function fixed_pts = FHNnc_fixed_points_multi(steady_state, options)
%Multiple fixed points
%   To obtain all the fixed points, look at all places where the 
%   V-nullcline crosses the w-nullcline

iter = -3:1e-1:3;
fixed_pts_temp = nan(size(iter));

for i = 1:length(iter)
    x1 = iter(i);
    [fixed_pts_temp(i), ~, exitflag_fixed_pts(i)] = fsolve(steady_state, x1, options);
    fixed_pts_temp_converged = exitflag_fixed_pts > 0;
    fixed_pts_temp_converged = fixed_pts_temp(fixed_pts_temp_converged);
end
fixed_pts = uniquetol(fixed_pts_temp_converged, 1E-3);
fixed_pts = sort(fixed_pts);

end

