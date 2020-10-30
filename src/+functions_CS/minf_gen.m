function [minf] = minf_gen(V)
    Vm_half = -29.735086419884503215196647032516;
    km = 9.4670757123944840690980421705083;

    minf = 1/(1 + exp((Vm_half - V)/ km));

end

