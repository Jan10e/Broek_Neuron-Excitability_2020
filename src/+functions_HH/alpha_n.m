function [alpha] = alpha_n(V)
    alpha =  0.01* ((10-V)./(exp((10-V)./10) -1));
end

