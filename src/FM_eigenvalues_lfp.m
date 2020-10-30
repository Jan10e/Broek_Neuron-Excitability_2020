% Title: Stability analysis of FHNm model - investigating Type I and II
% excitability. Get bifurcation points mathematically: where the tangent of
% the nullclines are equal
%
% Author: Jantine Broek
% Created: 21 Nov 2018
%%
clear
clc
dbstop if error

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020';

dir_work = '/src';
dir_data = '/data';
dir_fig = '/figures';

cd(fullfile(dir_base, dir_work));

%% functions_FM and Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));
w_inf3 = @(v, v0, w0) 2./(1 + exp(-5 * (v - v0))) + w0;

%parameters
v0 = 0.81; 
w0 = 0.01;

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');


%% (1.) Get SN coordinates
% det(jacobian) = 0
[v_sn, w_sn, I_app_sn] = functions_FM.FM_sn_coordinates(v0, w0, options);

I_app = I_app_sn - eps;

%% (2.) Get left fixed point just before SN bifurcation 
steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
    - ( w_inf1(y(1) - v0) + w0 ).^2 + I_app);

fixed_pts = functions_FM.FM_fixed_points_multi(steady_state, options);

v_lfp = min(fixed_pts);
w_lfp = w_inf3(v_lfp,v0,w0);

%% (3.) Eigenvalues for lfp for different epsilon
epsilon = 0:1e-3:0.5;

real_lambda=zeros(2,length(epsilon));
imag_lambda=zeros(2,length(epsilon));
for i = 1:length(epsilon)
    
    jacobian = [ [1 - v_lfp^2, -2 * w_lfp]; ...
        [epsilon(i) * (10 * exp(5*v0 - 5*v_lfp) ...
        / (exp(5*v0 - 5*v_lfp) + 1)^2), -epsilon(i)] ];
    
    real_lambda(:,i) = real(eig(jacobian));
    imag_lambda(:,i) = imag(eig(jacobian));
end


%% (4.) Plot real and imaginary values
f1 = figure;

subplot(2,1,1)
plot(epsilon, real_lambda,  'linewidth', 5)
ylim([-0.3 0.3])
xlabel('\epsilon'); 
ylabel('Re(\lambda)')
title(['Eigenvalues', '\newline', ...
    'V_0 = ', num2str(v0), ', w_0 = ', ...
    num2str(w0)], 'fontsize', 12);
hold on

plot(epsilon,zeros(length(epsilon)),'--k')
hold off


subplot(2,1,2)
plot(epsilon, imag_lambda, 'LineWidth', 5)
xlabel('\epsilon'); ylabel('Im(\lambda)')
hold on

plot(epsilon,zeros(length(epsilon)),'--k')


%% Export/Save
outfile = ['paper_fig4_FM_eigenvalues_lfp_example_V', num2str(v0), '_w', num2str(w0)];

suffix_fig = '';

out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% figures
%saveas(f1, outpath_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')

