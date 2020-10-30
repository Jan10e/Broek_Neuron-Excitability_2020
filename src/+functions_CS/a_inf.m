function [a] = a_inf(V)
% (Dayan-Abbott, p. 224)

a = (0.0761*exp(0.0314*(V+94.22))./(1+exp(0.0346*(V+1.17)))).^(1/3.0);

end

