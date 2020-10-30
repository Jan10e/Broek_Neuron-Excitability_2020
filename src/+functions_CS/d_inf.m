function [d,tau_d] = d_inf(V)

d = 1./(exp(-0.15.*(V+50))+1);
tau_d = 235./100; %*10 for ultraslow, *0.01 for fast