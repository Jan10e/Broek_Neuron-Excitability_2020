function [ht,tau_ht] = ht_inf(V)

ht = 1./(1+exp((V+31)./4.03)); %0.403 instead of 4.03, 51 instead of 81 for IVcurve, 31 for impulse

if V < -80
    tau_ht = 10*exp((V+467)./66.6); %none
else
    tau_ht = 10*exp(-(V+21.88)./10.2)+28;
end