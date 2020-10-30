function [mAinf] = mAinf_gen(V)
    VmA_half = -76.117031458553520175549113545954;
    kmA = 51.730805414763748477626424317013;

    mAinf = 1/(1 + exp((VmA_half - V)/ kmA));

end

