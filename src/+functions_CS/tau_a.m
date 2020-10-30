function [tau] = tau_a(V)
% (Dayan-Abbott, p. 224)

tau = 0.3632 + 1.158./(1+exp(0.0497*(V+55.96)));

end

