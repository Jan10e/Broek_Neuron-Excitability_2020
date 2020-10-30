function [ninf60] = ninf60_genV0(v, Vn_half)
    kn = 16.345749588434141227890436963931;
   
    ninf60 = 1./(1 + exp((Vn_half - v)./ kn));   
end

