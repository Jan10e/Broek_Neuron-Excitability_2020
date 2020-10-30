function h = hinf_ab(v)

h = functions.alpha_h(v)/(functions.alpha_h(v)+functions.beta_h(v));

% h = (7.*exp(- V./20 - 12./5))./(100.*((7.*exp(- V./20 - 12./5))./100 + 1./(exp(- V/10 - 9/5) + 1)));
end