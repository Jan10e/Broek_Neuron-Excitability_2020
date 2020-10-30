function [taum60] = taum60_ab(V)
    taum60 = 1/(functions.alpha_m(V+60) + functions.beta_m(V+60));
end


