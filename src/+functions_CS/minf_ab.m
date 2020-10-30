function m = minf_ab(v)

m = functions.alpha_m(v)/(functions.alpha_m(v)+functions.beta_m(v));

% m = -(V./10 + 297./100)./((exp(- V./10 - 297./100) - 1).*(4.*exp(- V./18 - 547./180) - (V./10 + 297./100)./(exp(- V./10 - 297./100) - 1)));
end