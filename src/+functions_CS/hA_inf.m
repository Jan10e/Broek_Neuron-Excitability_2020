function [hA] = hA_inf(V)

hA = 1./(1+exp((V+53.3)./14.54)).^4;

end

