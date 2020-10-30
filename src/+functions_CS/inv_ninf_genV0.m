function [inv_ninf] = inv_ninf_genV0(V, Vn_half)
    
    inv_ninf = Vn_half - (4600919484722715*log(1/V - 1))/281474976710656;
    
end

