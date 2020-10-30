function d = taud60(V)
%value for tau_d
d = 17 * (exp(-(V + 45)^2/600)) + 1.5;

end
