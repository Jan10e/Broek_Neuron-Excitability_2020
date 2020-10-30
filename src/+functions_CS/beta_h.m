function [beta] =beta_h(V)
    beta = 3.8/(1+exp(-0.1*(V+18)));
end

