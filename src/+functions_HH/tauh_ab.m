function [tauh] = tauh_ab(V)
    tauh = 1/(functions.alpha_h(V) + functions.beta_h(V));
end


