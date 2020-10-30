function [taum] = taum_ab(V)
    taum = 1/(functions.alpha_m(V) + functions.beta_m(V));
end


