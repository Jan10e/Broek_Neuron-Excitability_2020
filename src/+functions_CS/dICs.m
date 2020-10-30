function [gf, gs, Istatic] = dICs(V, gA, gCa)

 Iapp = -10;
   
%Reversal potentials
    VNa = 55;
    VK = -75;
    VCa = 120;
    Vl = -17;

    % conductances
    gl = 0.3;
    gNa = 120;
    gKDR = 20;
    taumCa = 235/100;

  gf = zeros(length(V));
  gs = zeros(length(V));
  Istatic = zeros(length(V));
  
  for i=1:length(V)
   wfs_m = var_contribution(taum(V(i)),taum(V(i)),taun(V(i)));
   wfs_h = var_contribution(tauh(V(i)),taum(V(i)),taun(V(i)));
   wfs_n = var_contribution(taun(V(i)),taum(V(i)),taun(V(i)));
   wfs_mA = var_contribution(taumA(V(i)),taum(V(i)),taun(V(i)));
   wfs_hA = var_contribution(tauhA(V(i)),taum(V(i)),taun(V(i)));
   wfs_mCa = var_contribution(taumCa,taum(V(i)),taun(V(i))); 
   
   gf(i) = -wfs_m*3*gNa*minf(V(i))^2*hinf(V(i))*(V(i)-VNa)*dotm(V(i)) -wfs_h*gNa*minf(V(i))^3*(V(i)-VNa)*doth(V(i)) - wfs_n*4*gKDR*ninf(V(i))^3*(V(i)-VK)*dotn(V(i)) - wfs_mA*3*gA*mAinf(V(i))^2*hAinf(V(i))*(V(i)-VK)*dotmA(V(i)) - wfs_hA*gA*mAinf(V(i))^3*(V(i)-VK)*dothA(V(i)) - 2*wfs_mCa*mCainf(V(i))*gCa*(V(i)-VCa)*dotmCa(V(i));
   gs(i) = -(1-wfs_m)*3*gNa*minf(V(i))^2*hinf(V(i))*(V(i)-VNa)*dotm(V(i)) -(1-wfs_h)*gNa*minf(V(i))^3*(V(i)-VNa)*doth(V(i)) -(1-wfs_n)*4*gKDR*ninf(V(i))^3*(V(i)-VK)*dotn(V(i)) - (1-wfs_mA)*3*gA*mAinf(V(i))^2*hAinf(V(i))*(V(i)-VK)*dotmA(V(i)) - (1-wfs_hA)*gA*mAinf(V(i))^3*(V(i)-VK)*dothA(V(i)) - 2*(1-wfs_mCa)*mCainf(V(i))*gCa*(V(i)-VCa)*dotmCa(V(i));

   Istatic(i) = -(-gNa*minf(V(i))^3*hinf(V(i))*(V(i)-VNa) -gKDR*ninf(V(i))^4*(V(i)-VK) - gA*mAinf(V(i))^3*hAinf(V(i))*(V(i)-VK) - gCa*mCainf(V(i))^2*(V(i)-VCa) - gl*(V(i)-Vl) + Iapp);

  end
  
end


