% Title: effect of epsilon on stability of left fixed point 
% Look at the eigenvalues of the lfp just before SN bifurcation
%
% Author: Jantine Broek
% Created: Oct 2019
%%
clear
clc

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/03_excitability/FHNnc/';

dir_work = '01_code/matlab';
dir_data = '02_data/data_matlab';
dir_fig = '03_figures/figures_matlab';

cd(fullfile(dir_base, dir_work));

%% Functions and Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));

v0 = 0.82; 

options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% (1.) Get SN coordinates
[v_sn, w_sn, I_app_sn] = functions.FHNnc_sn_coordinates(v0, options);

I_app = I_app_sn - eps;

%% (2.) Get left fixed point just before SN
steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
    - ( w_inf1(y(1)-v0)) + I_app);

fixed_pts = functions.FHNnc_fixed_points_multi(steady_state, options);
v_lfp = min(fixed_pts);

%% (3.) Real and Im values of lfp for different epsilon values
epsilon = 0:1e-3:0.5;

real_lambda=zeros(2,length(epsilon));
imag_lambda=zeros(2,length(epsilon));
for i = 1:length(epsilon)
    
    jacobian = [ [1-v_lfp^2, -1]; ...
        [(10*epsilon(i)*exp(5*v0 - 5*v_lfp))...
        /(exp(5*v0 - 5*v_lfp) + 1)^2, -epsilon(i)] ];
    
    lambda = eig(jacobian);
    
    real_lambda(:,i) = real(eig(jacobian));
    imag_lambda(:,i) = imag(eig(jacobian));
end

%% (4.) Plot real and imaginary values
f1 = figure;

subplot(2,1,1)
plot(epsilon, real_lambda,  'linewidth', 10)
xlabel('\epsilon'); 
ylabel('Re(\lambda)')
title(['Eigenvalues', '\newline', ...
    'V_0 = ', num2str(v0)], 'fontsize', 12);
hold on

plot(epsilon,zeros(length(epsilon)),'--k')
hold off


subplot(2,1,2)
plot(epsilon, imag_lambda, 'LineWidth', 10)
xlabel('\epsilon'); ylabel('Im(\lambda)')
hold on

plot(epsilon,zeros(length(epsilon)),'--k')

%% Export/Save
outfile = ['paper_fig3_FHNnc_eigenvalues_lfp_exmaple_V', num2str(v0)];

suffix_fig = '';

% out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

% outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% Figures
% saveas(f1, outpath_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')
