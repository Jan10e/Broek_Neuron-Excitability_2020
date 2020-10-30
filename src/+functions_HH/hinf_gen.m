function [hinf, hinf60] = hinf_gen(V)
    Vh_half = 2.6924982163435065908509849483864;
    kh = -6.9501544723237258866654620505739;
   
    hinf = 1/(1 + exp((Vh_half - V)/ kh));   
end


