function [m] = dotm(V)
m = - 1/(10*(exp(- V/10 - 297/100) - 1)*(4*exp(- V/18 - 547/180) - (V/10 + 297/100)/(exp(- V/10 - 297/100) - 1))) - (exp(- V/10 - 297/100)*(V/10 + 297/100))/(10*(exp(- V/10 - 297/100) - 1)^2*(4*exp(- V/18 - 547/180) - (V/10 + 297/100)/(exp(- V/10 - 297/100) - 1))) - ((V/10 + 297/100)*((2*exp(- V/18 - 547/180))/9 + 1/(10*(exp(- V/10 - 297/100) - 1)) + (exp(- V/10 - 297/100)*(V/10 + 297/100))/(10*(exp(- V/10 - 297/100) - 1)^2)))/((exp(- V/10 - 297/100) - 1)*(4*exp(- V/18 - 547/180) - (V/10 + 297/100)/(exp(- V/10 - 297/100) - 1))^2);
end

