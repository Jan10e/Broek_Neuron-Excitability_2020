function [hinf60] = hinf60_ab(v)
    hinf60 = functions.alpha_h(v + 60) ./(functions.alpha_h(v + 60) + functions.beta_h(v + 60));
end

