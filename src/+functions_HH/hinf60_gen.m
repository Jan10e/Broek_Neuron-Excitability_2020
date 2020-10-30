function [hinf60] = hinf60_gen(V)
    Vh_half = -57.307501783656493409149015051614;
    kh = -6.9501544723237258866654620505739;
   
%     hinf60 = 1/(1 + exp((Vh_half - (V+60))/ kh));   
    hinf60 = 1/(1 + exp((Vh_half - V)/ kh));   
end


