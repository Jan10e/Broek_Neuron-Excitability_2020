function [v, w] = FHNnc_forward_euler(v0, epsilon, I_app, v_init, w_init, t)
% Simulation of mirrored Fitzhugh Nagumo model

dt = t(2)-t(1);

%% Functions
w_inf = @(v) 2./(1 + exp(-5 * (v)));

%% Euler solve info
v = v_init;
w = w_init;

%% Integrate using Forward Euler
    for i = 1:length(t)-1
         w(i + 1) = w(i) + dt .* (epsilon*(w_inf(v(i) - v0) - w(i)));
         v(i + 1) = v(i) + dt .* (v(i) - v(i)^3/3 - w(i) + I_app);
    end

end

