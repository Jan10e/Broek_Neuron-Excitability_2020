function [minf60] = minf60_genV0(V, Vm_half)
    km = 9.4718463719810198256478298087987;
   
%     minf60 = 1/(1 + exp((Vm_half - (V + 60))./ km));   
    minf60 = 1/(1 + exp((Vm_half - V)./ km));   
end

