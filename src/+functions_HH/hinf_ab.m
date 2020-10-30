function [hinf] = hinf_ab(V)
    hinf = functions.alpha_h(V) ./(functions.alpha_h(V) + functions.beta_h(V));
end

