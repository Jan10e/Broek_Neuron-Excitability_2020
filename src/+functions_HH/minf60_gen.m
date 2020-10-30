function [minf60] = minf60_gen(v)
    Vm_half = -35.024577352365961579525247935964;
    km = 9.4718463719810198256478298087987;
   
%     minf60 = 1/(1 + exp((Vm_half - (V + 60))./ km));   
    minf60 = 1/(1 + exp((Vm_half - v)./ km));   
end



