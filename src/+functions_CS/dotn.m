function [n] = dotn(V)
n = - 1/(100*(exp(- V/10 - 457/100) - 1)*(exp(- V/80 - 557/800)/8 - (V/100 + 457/1000)/(exp(- V/10 - 457/100) - 1))) - (exp(- V/10 - 457/100)*(V/100 + 457/1000))/(10*(exp(- V/10 - 457/100) - 1)^2*(exp(- V/80 - 557/800)/8 - (V/100 + 457/1000)/(exp(- V/10 - 457/100) - 1))) - ((V/100 + 457/1000)*(exp(- V/80 - 557/800)/640 + 1/(100*(exp(- V/10 - 457/100) - 1)) + (exp(- V/10 - 457/100)*(V/100 + 457/1000))/(10*(exp(- V/10 - 457/100) - 1)^2)))/((exp(- V/10 - 457/100) - 1)*(exp(- V/80 - 557/800)/8 - (V/100 + 457/1000)/(exp(- V/10 - 457/100) - 1))^2);
end

