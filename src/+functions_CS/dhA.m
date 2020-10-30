function hA = dhA(V, hA)

hA = (1/functions.tau_hA(V))*(functions.hA_inf(V) - hA);

end

