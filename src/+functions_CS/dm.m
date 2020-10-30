function [m] = dm(v, m)

m =(1/functions.tau_m(v))*(functions.minf_ab(v) - m);

end

