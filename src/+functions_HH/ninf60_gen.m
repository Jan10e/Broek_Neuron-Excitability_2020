function [ninf60] = ninf60_gen(V)
    Vn_half = -48.412914431643244507103066115362;
    kn = 16.345749588434141227890436963931;
   
%     ninf60 = 1/(1 + exp((Vn_half - (V + 60))./ kn));   
    ninf60 = 1./(1 + exp((Vn_half - V)./ kn));   
end

