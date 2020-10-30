function [tau] = tau_n(V)

tau = 1/(functions.alpha_n(V) + functions.beta_n(V)); 

% tau = 2./((19.*exp(- V./80 - 557/800))./40 - (19.*(V./100 + 457./1000))/(5.*(exp(- V./10 - 457./100) - 1)));
end

