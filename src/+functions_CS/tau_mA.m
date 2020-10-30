function [tau] = tau_mA(V)
tau = 1*(0.3632 + 1.158./(1+exp((V+55.96)./20.12)));
end

