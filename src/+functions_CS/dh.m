function h = dh(v, h)

h = (1/functions.tau_h(v))*(functions.hinf_ab(v) - h);

end

