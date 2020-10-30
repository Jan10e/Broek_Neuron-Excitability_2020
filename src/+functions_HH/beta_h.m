function [beta] = beta_h(V)
    beta = 1/(exp((30-V)/10) + 1);
end

