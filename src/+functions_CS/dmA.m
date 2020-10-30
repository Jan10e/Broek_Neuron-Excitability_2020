function mA = dmA(v, mA)

mA= (1/functions.tau_mA(v))*(functions.mA_inf(v) - mA);

end

