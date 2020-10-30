function [b] = b_inf(V)
% (Dayan-Abbott, p. 224)

b = (1./(1+exp(0.0688*(V+53.3)))).^4;

end

