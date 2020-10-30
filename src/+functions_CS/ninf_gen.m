function [ninf] = ninf_gen (V)
    Vn_half = -44.112914431643244507103066115362;
    kn = 16.345749588434141227890436963931;


    ninf = 1./(1 + exp((Vn_half - V)./ kn));
    
end

