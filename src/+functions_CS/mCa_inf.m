function [mCa] = mCa_inf(V)
mCa = 1./(exp(-0.15.*(V+50))+1);
end

