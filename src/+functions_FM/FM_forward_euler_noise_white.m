function y = FM_forward_euler_noise_white(V0, w0, epsilon, I_app, V_init, w_init, t, sigma)
%simulation of mirrored Fitzhugh Nagumo model

dt = t(2)-t(1);

%% Functions
w_inf1 = @(V) 2./(1 + exp(-5 * (V)));

y = zeros(2, length(t));
y(:,1) = [V_init; w_init];

%% Generate white noise Vector
eta = sigma.*randn(length(t),1).*sqrt(dt); 

%% Integrate using Forward Euler
for i=1:numel(t)-1
    F = @(t,y) [ (y(1) - y(1).^3/3 - y(2)^2 + I_app) +eta(i); epsilon*(w_inf1(y(1) - V0) + (w0 - y(2))) ];
    y(:,i + 1) = y(:,i) + dt* F(t, y(:,i));
  
end

end

