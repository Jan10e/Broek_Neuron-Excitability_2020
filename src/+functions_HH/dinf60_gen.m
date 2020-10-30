function [dinf] = dinf60_gen(V)
    Vd_half = -55.0;
    kd = 3.0;

    dinf = 1/(1 + exp((Vd_half - V)/ kd));
end

