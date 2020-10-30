function [ninf] = ninf_gen(V)
    Vn_half = 11.587085568356755492896933884638;
    kn = 16.345749588434141227890436963931;
   
    ninf = 1/(1 + exp((Vn_half - V)/ kn));   
end

