function [m,tau_m] = mKslow_inf(V)

NSHFT=-14.3;

alpha_m = -0.01.*(V+50+NSHFT)./(exp(-(V+50+NSHFT)./10)-1);

beta_m = 0.125.*exp(-(V+60+NSHFT)./80);


a = alpha_m + beta_m;

m = alpha_m ./ a;
tau_m =100 ./ (3.8.*a);