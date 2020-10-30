function [hAinf] = hAinf_gen(V)
    VhA_half = -77.507835272251688786779672864736;
    khA = -11.423375550577798153874266558853;

    hAinf = 1/(1 + exp((VhA_half - V)/ khA));

end

