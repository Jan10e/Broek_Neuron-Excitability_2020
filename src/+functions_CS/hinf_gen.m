function [hinf] = hinf_gen(V)
    Vh_half = -45.307501783656493409149015051614;
    kh = -6.9501544723237258866654620505739;
    
    hinf = 1/(1 + exp((Vh_half - V)/ kh));

end

