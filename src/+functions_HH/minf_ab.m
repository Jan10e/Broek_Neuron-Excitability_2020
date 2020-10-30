function [minf] = minf_ab(V)
    minf = functions.alpha_m(V)./(functions.alpha_m(V) + functions.beta_m(V));
end

