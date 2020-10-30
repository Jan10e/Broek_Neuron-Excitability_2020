function [ainf] = ainf_gen(V)
    Va_half = -76.118726360804529569454690710534;
    ka = 51.737193395598845206959653408111;
    
    ainf = 1/(1 + exp((Va_half - V)/ ka));

end

