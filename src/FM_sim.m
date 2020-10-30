% Title: Simulation of mirrored FHN model - investigating Type I and II
% excitability
%
% Author: Jantine Broek
% Created: 26 Nov 2018
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

v0 = 0.81;
w0 = 0.01; 
epsilon = 2e-1;
% I_app = 0.66677;    % use this for fixed I_app; comment the 'I_SN: for situations right before or after SN'

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');

%% I_SN: for situations right before or after SN
[v_sn, w_sn, I_app_sn] = functions_FM.FM_sn_coordinates(v0, w0, options);

if ~isnan(w_sn)
    I_app = I_app_sn - 1e-5;
%     I_app = I_app_sn;
else
    I_app = 0.81;
end

%% Forward Euler solve info

% integration time
tend = (5/epsilon)*10; %make sure to capture some spikes
% tend = 1e4; %2e3 for 0.5Hz; 1e4 for 0.1Hz 
dt = 0.01;
t = 0:dt:tend;
% Ttrans = 0.3/dt;
% Ttrans_idx = find(t == Ttrans); 
Hz = 1000/tend;

% initial values
v_init = -1.2; 
% v_init = 0; %to look for bistability
w_init = w0; %generally start at w0


%% Integrate: forward Euler
[v, w] = functions_FM.FM_forward_euler(v0, w0, epsilon, ...
    I_app, v_init, w_init, t);

%% Plot V(m) vs time
% f1 = figure('Position',[800 200 700 700]);
f1 = figure('Position',[1000 200 500 700]);
subplot(2, 1, 1);

plot(t, v, 'k-', 'linewidth', 4);
% plot(t(Ttrans_idx:end), V(Ttrans_idx:end), 'k', 'LineWidth', 1.5)
ylim([-3 3])
xlabel('$time (ms)$', 'interpreter','latex','fontsize',11); 
ylabel('$V$', 'interpreter','latex','fontsize',11);
title(['Simulation', '\newline', ...
    'V_0 = ', num2str(v0), ', w_0 = ', num2str(w0), ...
    ', \epsilon = ', num2str(epsilon), '\newline', ...
    'I_{app} = ', num2str(I_app), ', freq = ', num2str(Hz), ' Hz'], ...
    'fontsize', 12);

%% Plot phase plane
subplot(2, 1, 2)
% f1 = figure;

% [vq,wq] = meshgrid(-3:.2:3, -3:.2:3); % quiver plot
% dv = vq - (1/3) * vq.^3 - wq.^2 + I_app;
% dw = epsilon * (w_inf1(vq - v0) + (w0 - wq));
% quiver(vq,wq,dv,dw, 1.5, 'Linewidth', 2)


% nullclines
syms vv ww 
dot_v = vv - (1/3) * vv^3 - ww^2 + I_app;
dot_w = w_inf1(vv - v0) + (w0 - ww);

hold on
nullcline_v = ezplot(dot_v, [-3 3]);
set(nullcline_v,'color',[0, 0.4470, 0.748], 'Linewidth', 5);
nullcline_w = ezplot(dot_w, [-3 3]);
set(nullcline_w, 'color', [0.4660, 0.6740, 0.1880], 'Linewidth', 5);

% trajectory
plot(v, w, 'color',[1 0.08 0.57],'LineWidth', 5) 

xlabel('$V$', 'interpreter','latex','fontsize',11);
ylabel('$w$', 'interpreter','latex','fontsize',11);
title(['Phase Plane', '\newline', ...
    'V_0 = ', num2str(v0), ', w_0 = ', num2str(w0), ...
    ', \epsilon = ', num2str(epsilon), ', I_{app} = ', ...
    num2str(I_app)], 'fontsize', 12);


%% Fixed points
% When three fixed point are present: SN
% When one fixed point is present: Hopf only

steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
    - ( w_inf1(y(1) - v0) + w0 ).^2 + I_app);


fixed_pts = functions_FM.FM_fixed_points_multi(steady_state, options);

%classication of fixed points
for c = 1:length(fixed_pts) 
    
    v_steady_state = fixed_pts(1,c); 
    w_steady_state = w_inf3(v_steady_state, v0, w0);
    
    jacobian =[ [1 - v_steady_state^2, -2 * w_steady_state]; ...
        [(epsilon*10*exp(5*v0 - 5*v_steady_state))...
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


[freq, ~, ~] = functions_FM.FM_freq(v0, w0, epsilon, ...
    I_app, v_init, w_init, t);

%% Export/Save
outfile = ['FM_sim_V', num2str(v0), '_w', num2str(w0), '_e',...
           num2str(epsilon), '_I', num2str(I_app)];
       
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
saveas(f1, outpath_fig_png,'png')

