function [minf] = minf_gen(V)
    Vm_half = 24.975422647634038420474752064036;
    km = 9.4718463719810198256478298087987;
   
    minf = 1/(1 + exp((Vm_half - (V))./ km));   
end



