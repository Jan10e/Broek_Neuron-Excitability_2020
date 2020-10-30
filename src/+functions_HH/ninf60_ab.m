function [ninf60] = ninf60_ab(V)
    ninf60 = functions.alpha_n(V + 60) ./(functions.alpha_n(V + 60) + functions.beta_n(V + 60));
end

