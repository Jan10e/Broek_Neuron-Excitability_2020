function [inv] = inv_ninf(V)
inv =  -(200*log(1./V - 21/20))/13 - 208/5;
end

