function [mt,tau_mt] = mt_inf(V)

mt = 1./(1+exp(-(V+47)./6.2)); %0.62 instead of 6.2, 47 instead of 67 for IVcurve 
tau_mt = 0.5*(0.612 + 1./(exp(-(V+131.6)./16.7)+exp((V+16.8)./18.2))); %0.5 added