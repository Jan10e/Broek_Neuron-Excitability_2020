function [alpha] = alpha_m(V)
    alpha = 0.1*((25-V)/(exp((25-V)/10) -1) );
end

