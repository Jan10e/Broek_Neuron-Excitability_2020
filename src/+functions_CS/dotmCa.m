function [mCa] = dotmCa(V)
mCa = (3*exp(- (3*V)/20 - 15/2))/(20*(exp(- (3*V)/20 - 15/2) + 1)^2);
end

