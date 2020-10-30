function [alpha] = alpha_n(V)
    alpha = 0.02*(V+45.7)/(1-exp(-0.1*(V+45.7)));
end

