function [tau] = tau_hA(V)
tau = 1*(1.24+2.678./(1+exp((V+50)./16.027)));
end

