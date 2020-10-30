function [v, w] = FM_forward_euler_winf_new2(v0, w0, epsilon, I_app, v_init, w_init, t)
%simulation of mirrored Fitzhugh Nagumo model

dt = t(2)-t(1);

%% Functions
w_inf1 = @(v) 10./(1 + exp(-5 * (v - (1/3))));

%% Euler solve info
v = v_init;
w = w_init;

%% Integrate using Forward Euler
    for i = 1:length(t)-1
         v(i + 1) = v(i) + dt .* (v(i) - v(i)^3/3 - w(i)^2 + I_app);
         w(i + 1) = w(i) + dt .* (epsilon*(w_inf1(v(i) - v0) + (w0 - w(i))));
    end

end
