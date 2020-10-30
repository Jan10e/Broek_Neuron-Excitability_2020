function [n] = ninf_ab(v)

n = functions.alpha_n(v)./(functions.alpha_n(v)+functions.beta_n(v)); 

% n = 1./(1.05+exp(0.065.*(-41.6-V)));

%-(V/50 + 457/500)/((exp(- V/10 - 457/100) - 1)*(exp(- V/80 - 557/800)/4 - (V/50 + 457/500)/(exp(- V/10 - 457/100) - 1)))
end