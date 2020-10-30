function [alpha] = alpha_h(V)
    alpha = 0.266*exp(-0.05*(V+48));
end

