function [ninf] = ninf_genV0(V, Vn_half)
    kn = 16.345749588434141227890436963931;

    ninf = 1./(1 + exp((Vn_half - V)./ kn));
    
end

