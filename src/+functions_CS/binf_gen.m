function [binf] = binf_gen(V)
    Vb_half = -77.4993171126280437153918549318;
    kb = -11.419355937287872822640697033497;

    binf = 1/(1 + exp((Vb_half - V)/ kb));
    
end

