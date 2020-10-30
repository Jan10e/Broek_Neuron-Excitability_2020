function [alpha] = alpha_m(V)
    alpha = 0.38*(V+29.7)/(1-exp(-0.1*(V+29.7)));
end

