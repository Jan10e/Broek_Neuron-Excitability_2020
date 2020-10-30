function tau = tau_m(v)

tau = 1/(functions.alpha_m(v) + functions.beta_m(v)); 

% tau = 1/((76*exp(- V/18 - 547/180))/5 - (19*(V/10 + 297/100))/(5*(exp(- V/10 - 297/100) - 1)));
end

