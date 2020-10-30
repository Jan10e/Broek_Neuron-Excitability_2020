function [inv_ninf] = inv_ninf_gen(V)
% this is the inverse function of the general eqn for ninf, i.e. 1/(1+ exp((Vhalf-V)/k))
    
 inv_ninf = -(4600919484722715*log(1/V - 1))/281474976710656 - 1552085195285743/35184372088832;
end

