% Title: Simulation of FHN model with a sigmoidal w-nullcline that 
% is the activation function.
%
% Author: Jantine Broek
% Created: Oct 2019
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

%% functions_FHNnc and Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));

v0 = 0.82; %canonical 0.07
epsilon = 0.004; %canonical 0.08
% I_app = 0.66789; %0.3 - 1.2; use this for a fixed I_app, but then comment
% out the `Get I_SN value' part

options = optimoptions('fsolve','Display','none', ...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14, ...
    'TolX',1e-14);

%% Get I_SN value
[v_sn, w_sn, I_sn] = functions_FHNnc.FHNnc_sn_coordinates(v0, options);

if ~isnan(w_sn)
    I_app = I_sn - 1e-5;
else 
    I_app = 0.71;
end

%% Integrate: Euler solve
tend = (5/epsilon)*10; %make sure to capture some spikes
% tend = 2e3; %2e3 for 0.5Hz; 
dt = 0.01;
t = 0:dt:tend;

v_init = -1.5; 
w_init = 0;

[v, w] = functions_FHNnc.FHNnc_forward_euler(v0, epsilon, I_app, v_init, w_init, t);

%% Plot V(m) vs time
f1 = figure('Position',[1000 200 500 700]);
subplot(2, 1, 1)

plot(t, v, 'k-', 'linewidth', 4);
xlabel('$time (ms)$', 'interpreter','latex','fontsize',12); 
ylabel('$V$', 'interpreter','latex','fontsize',12);
ylim([-3 3])
title(['V0 = ', num2str(v0), ', \epsilon = ', ...
    num2str(epsilon), ', I_{app} = ', num2str(I_app)], ...
    'fontsize', 12);

%% Plot phase plane
subplot(2, 1, 2)
% figure;

% [vq,wq] = meshgrid(-3:0.4:3, -3:0.4:3); % quiver plot
% dv = vq - (1/3) .* vq.^3 - wq + I_app;
% dw = epsilon * (w_inf1(vq - v0) - wq);
% quiver(vq, wq, dv, dw, ...
%      1.5,'color', 'black', 'Linewidth', 1) %1.5 is the scale factor

% nullclines
syms vv ww 
dot_v = vv - (1/3) * vv^3 - ww + I_app;
dot_w = epsilon * (w_inf1(vv - v0) - ww);

hold on
nullcline_v = ezplot(dot_v, [-3 3]);
set(nullcline_v,'color',[0, 0.4470, 0.7410],'LineWidth',5);
nullcline_w = ezplot(dot_w, [-3 3]);
set(nullcline_w,'color',[0.4660, 0.6740, 0.1880],'LineWidth',5);

% trajectory
plot(v, w,'color',[1 0.08 0.57],'LineWidth', 5); 

xlabel('$V$', 'interpreter','latex','fontsize',12);
ylabel('$w$', 'interpreter','latex','fontsize',12);
title(['Phase Plane', '\newline', ...
    'V_0 = ', num2str(v0), ', \epsilon = ', ...
    num2str(epsilon), ', I_{app} = ', num2str(I_app)], ...
    'fontsize', 12);


%% Fixed points
% When three fixed point are present: SN
% When one fixed point is present: Hopf only

steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
    - ( w_inf1(y(1)-v0)) + I_app);

fixed_pts = functions_FHNnc.FHNnc_fixed_points_multi(steady_state, options);

for c = 1:length(fixed_pts) %classication of fixed points
    
    v_steady_state = fixed_pts(1,c); 
    w_steady_state = w_inf1(v_steady_state - v0);
    
    jacobian = [ [1 - v_steady_state^2, -1]; ...
        [(10*epsilon*exp(5*v0 - 5*v_steady_state))...
        /(exp(5*v0 - 5*v_steady_state) + 1)^2, -epsilon] ];
    
    lambda(:,c) = eig(jacobian);

    if imag(lambda(1,c))~= 0 && real(lambda(1,c))<0 ...
            && real(lambda(2,c))<0
        fprintf('fixed point %d, %d is stable spiral.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15, 'LineWidth', 3)
    end
    if imag(lambda(1,c))~=0 && real(lambda(1,c))>0 ...
            && real(lambda(2,c))>0
        fprintf('fixed point %d, %d is unstable spiral.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
    end
    if imag(lambda(1,c))==0 && real(lambda(1,c))<0 ...
            && real(lambda(2,c))<0
        fprintf('fixed point %d, %d is stable node.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15, 'LineWidth', 3)
    end
    if imag(lambda(1,c))==0 && real(lambda(1,c))>0 ...
            && real(lambda(2,c))>0
        fprintf('fixed point %d, %d is unstable node.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
    end
    if imag(lambda(1,c))==0 && real(lambda(1,c))>0 ...
            && real(lambda(2,c))<0
        fprintf('fixed point %d, %d is a saddle.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'x', 'color', [1.0 0.6 0.4], ...
            'MarkerSize', 40, 'LineWidth', 5)
    end
    if imag(lambda(1,c))==0 && real(lambda(1,c))<0 ...
            && real(lambda(2,c))>0
        fprintf('fixed point %d, %d is a saddle.\n', v_steady_state, w_steady_state);
        plot(v_steady_state, w_steady_state, 'x', 'color', [1.0 0.6 0.4], ... 
            'MarkerSize', 40, 'LineWidth', 5)
    end
end  


%% Export/Save
outfile = ['paper_fig3_FHNnc_sim_V', num2str(v0), '_e', num2str(epsilon), ...
    '_I', num2str(I_app)];

suffix_fig = '';
suffix_data = '';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);


% figures
% saveas(f1, outpath_fig_eps,'eps')
% saveas(gcf, outpath_fig_png,'png')
