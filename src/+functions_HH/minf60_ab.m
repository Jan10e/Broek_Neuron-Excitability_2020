function [minf60] = minf60_ab(v)
    minf60 = functions.alpha_m(v + 60)./(functions.alpha_m(v + 60) + functions.beta_m(v + 60));
end

