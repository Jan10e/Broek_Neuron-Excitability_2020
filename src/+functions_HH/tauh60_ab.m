function [tauh60] = tauh60_ab(v)
    tauh60 = 1/(functions.alpha_h(v + 60) + functions.beta_h(v + 60));
end


