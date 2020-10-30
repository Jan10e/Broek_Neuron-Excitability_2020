function mCa = dmCa(v, mCa)
tau_mCa = 235/100;

mCa =(1/tau_mCa)*(functions.mCa_inf(v) - mCa);

end

