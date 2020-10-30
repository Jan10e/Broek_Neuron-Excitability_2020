function [taun60] = taun60_ab(V)
    taun60 = 1/(functions.alpha_n(V + 60) + functions.beta_n(V + 60));
end

