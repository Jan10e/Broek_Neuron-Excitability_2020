function [ninf] = ninf_ab(v)
    ninf = functions.alpha_n(v)./(functions.alpha_n(v) + functions.beta_n(v));
end

