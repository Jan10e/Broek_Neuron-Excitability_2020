function tau = tau_h(v)

tau = 1/(functions.alpha_h(v) + functions.beta_h(v)); 

% tau = 1/((133*exp(- V/20 - 12/5))/500 + 19/(5*(exp(- V/10 - 9/5) + 1)));
end

