function [wfs] = var_contribution(tauX, tauf, taus)
%Computes contribution of time constant to either fast or slow timescale

wfs = 1;

  if tauX < tauf
   wfs = 1;
  elseif tauf <= tauX < taus
   wfs = (log(taus)-log(tauX))/(log(taus)-log(tauf));
  else
   wfs = 0;
  end

end
