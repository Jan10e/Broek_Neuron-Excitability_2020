function [dinf] = dinf_gen(V)
    Vd_half = 5.0;
    kd = 3.0;

    dinf = 1/(1 + exp((Vd_half - V)/ kd));
end


