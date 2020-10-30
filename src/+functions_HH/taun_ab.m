function [taun] = taun_ab(V)
    taun = 1/(functions.alpha_n(V) + functions.beta_n(V));
end

