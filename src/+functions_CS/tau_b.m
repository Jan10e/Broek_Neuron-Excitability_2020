function [tau] = tau_b(V)
% (Dayan-Abbott, p. 224)

tau = 1.24 + 2.678./(1+exp(0.0624*(V+50)));

end

