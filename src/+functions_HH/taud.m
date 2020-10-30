function d = taud(V)
%value for tau_d
d = 17 * (exp(-(V - 15)^2/600)) + 1.5;

end