function [m,tau_m] = mH_infN(V)

 m = 1./(1+exp((V+80)./6));
 tau_m = 1*(272 + 1499./(1+exp((V+42.2)./(-8.73))));

end

