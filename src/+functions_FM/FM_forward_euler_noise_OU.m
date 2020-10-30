function [v, w] = FM_forward_euler_noise_OU(v0, w0, epsilon, I_app, v_init, w_init, t, sigma, tau_noise)
%simulation of mirrored Fitzhugh Nagumo model
%some reference: http://www.scholarpedia.org/article/Stochastic_dynamical_systems

dt = t(2)-t(1);
OU = zeros(1,length(t));

%% Functions
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));

dv = @(v, w, I_app, t) v - v^3/3 - w^2 + I_app;
dw = @(v, w, v0, w0, epsilon, t) epsilon*(w_inf1(v - v0) + (w0 - w));


%% Euler solve info
v = v_init;
w = w_init;

%% Generate Ornstein-Uhlenbeck (OU) noise vector
for i = 1:length(t)-1
    OU(i+1) = OU(i) + (-OU(i) * dt + sqrt(2 .* dt * tau_noise * sigma.^2).*(randn)) *(1./tau_noise);
end


%% Integrate using Forward Euler
for ii = 1:length(OU)-1

    v(ii + 1) = v(ii) + dt * dv(v(ii), w(ii), I_app + OU(ii), t);
    w(ii + 1) = w(ii) + dt * dw(v(ii), w(ii), v0, w0, epsilon, t);
 
end

end

