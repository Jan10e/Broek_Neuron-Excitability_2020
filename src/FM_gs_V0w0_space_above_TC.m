% Title: Get the g_s(V_th) value at epsilon_crit. Separately get the dV/dw and
% dw_infty/dV, as gs = dV/dw * dw_infty/dV for all critical values above
% TC.
%
% Author: Jantine Broek
% Created: Dec 2019
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

%% Loop over parameter space
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));
w_inf3 = @(v, v0, w0) 2./(1 + exp(-5 * (v - v0))) + w0;

v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;
epsilon_range = 0:1e-5:0.4;

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');


%% Get TC coordinates
w0_TC = functions_FM.FM_TC_line(v0_range, options);

%% Specify an Iapp for the case of Hopf only
v0_boundary = -0.42;
w0_boundary = -0.05;

% Get SN coordinates at boundary
[v_sn_bdry, w_sn_bdry, I_app_sn_bdry] = functions_FM.FM_sn_coordinates(v0_boundary, w0_boundary, options);
I_app_bdry = I_app_sn_bdry - eps;


%% Bifurcation map loop
tic
for p = 1:length(v0_range)
    v0 = v0_range(p)
    
    for q = 1:length(w0_range)
        w0 = w0_range(q);
        
        %% (1.) Get SN coordinates
        [v_sn, w_sn, I_app_sn] = functions_FM.FM_sn_coordinates(v0, w0, options);
        
        if ~isnan(w_sn) %three fixed points present > SN
            
            %% (2.) Get left fixed points before SN
            I_app = I_app_sn - eps;
            
            steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
                - ( w_inf1(y(1) - v0) +w0 ).^2 + I_app);
            
            fixed_pts = functions_FM.FM_fixed_points_multi(steady_state, options);
            v_lfp = min(fixed_pts);
            w_lfp = w_inf3(v_lfp,v0,w0);
            
            
            %% (3.) Get real values of lfp
            for r = 1:length(epsilon_range)
                epsilon = epsilon_range(r);
                
                jacobian = [ [1 - v_lfp^2, -2 * w_lfp]; ...
                    [epsilon * (10 * exp(5*v0 - 5*v_lfp) ...
                    / (exp(5*v0 - 5*v_lfp) + 1)^2), -epsilon] ];
                
                real_lambda = real(eig(jacobian));
                
                
                if min(real_lambda) < 0 %select epsilon values of negative Re(lambda) values: SNIC
                    epsilon_crit(p,q) = epsilon_range(r);
                    pv_pw(p,q) = 2 * w_inf3(v_lfp,v0,w0);
                    pwinf_pv(p,q) = ( 10*exp(5*v0 - 5*v_lfp))...
                        /(exp(5*v0 - 5*v_lfp) + 1)^2;
                    gs_vth(p,q) = ( 2*(w_inf3(v_lfp, v0, w0)) )...
                        .* ( (10*exp(5*v0 - 5*v_lfp))/(exp(5*v0 ...
                        - 5*v_lfp) + 1)^2 );
                    break
                end
            end
            
            
        elseif isnan(w_sn) %one fixed point present
            
            %% (1.) Get fixed point
            steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
                - ( w_inf1(y(1) - v0) +w0 ).^2 + I_app_bdry);
            
            fixed_pt = functions_FM.FM_fixed_points_multi(steady_state, options);
            v_Hopf = min(fixed_pt);
            w_Hopf = w_inf3(v_Hopf,v0,w0);
            
            
            %% (2.) Get real values of fp
            for r = 1:length(epsilon_range)
                epsilon = epsilon_range(r);
                
                jacobian = [ [1 - v_Hopf^2, -2 * w_Hopf]; ...
                    [epsilon * (10 * exp(5*v0 - 5*v_Hopf) ...
                    / (exp(5*v0 - 5*v_Hopf) + 1)^2), -epsilon] ];
                
                real_lambda = real(eig(jacobian));
                
                
                if min(real_lambda) < 0 %select epsilon values of negative Re(lambda) values
                    epsilon_Hopf_crit(p,q) = epsilon_range(r);
                    pv_pw_Hopf(p,q) = 2 * w_inf3(v_Hopf,v0,w0);
                    pwinf_pv_Hopf(p,q) = ( 10*exp(5*v0 - 5*v_Hopf))...
                        /(exp(5*v0 - 5*v_Hopf) + 1)^2;
                    gs_vth_Hopf(p,q) = ( 2*(w_inf3(v_Hopf, v0, w0)) )...
                        .* ( (10*exp(5*v0 - 5*v_Hopf))/(exp(5*v0 ...
                        - 5*v_Hopf) + 1)^2 );
                    break
                end
            end
                
            
        end
    end
end
toc


%% Plot figures
%plot bifurcation map for teps_crit
f1 = figure;

subplot(4,1,1)
colormap(jet)
% imagesc(w0_range, v0_range, epsilon_crit)
imagesc(w0_range, v0_range, log10(epsilon_crit))    
hold on

% plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title('$\epsilon_{crit}$','interpreter','latex','fontsize',12)


%plot bifurcation map for gs
subplot(4,1,2)
colormap(jet)
imagesc(w0_range, v0_range, abs(log10(gs_vth)))
% imagesc(w0_range, v0_range, gs_crit_pos)
hold on

% plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title('$g_s(V_{th})$','interpreter','latex','fontsize',12)



%plot bifurcation map for pv_pw
subplot(4,1,3)
colormap(jet)
imagesc(w0_range, v0_range, abs(log10(pv_pw)))
% imagesc(w0_range, v0_range, pv_pw_pos)
hold on

% plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title('$\frac{\partial \dot{V}}{\partial w}(V_{th})$', 'interpreter','latex','fontsize',12)



%plot bifurcation map for pwinf_pv
subplot(4,1,4)
colormap(jet)
imagesc(w0_range, v0_range, log10(pwinf_pv))
% imagesc(w0_range, v0_range, pwinf_pv_pos)
hold on

% plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title('$\frac{\partial w_\infty}{\partial V}(V_{th})$', 'interpreter','latex','fontsize',12)


%% gs vs epsilon (this one takes a while)
% f2 = figure;
% hold on
% for i = 1:size(epsilon_crit, 1)
%     for j = 1:size(epsilon_crit, 2)
%         plot(gs_vth(i,j), epsilon_crit(i,j), 'ok')
%         hold on
%         plot([0, 0.3], [0, 0.3], 'g-', 'Linewidth', 1.5)
%     end
% end
% % xlim([0, 1.4])
% xlabel('$g_s(V_{th})$', 'interpreter','latex','fontsize',12);
% ylabel('$\epsilon$', 'interpreter','latex','fontsize',12);
% title('\epsilon_{crit} vs g_s(V_{th})')


%% Export/Save
outfile = 'FM_epsilon_crit_V0w0_space_HS_gs';

suffix_data_gs_pos = '_gs_vth_pos';
suffix_data_pvpw_pos = '_pvpw_pos';
suffix_data_pwinfpv_pos = '_pwinfpv_pos';
suffix_data_gs_original = '_gs_vth_original';
suffix_data_pvpw_original = '_pvpw_original';
suffix_data_pwinfpv_original = '_pwinfpv_original';
suffix_data_epsilon = '_epsilon_crit';

suffix_fig1 = '';
% suffix_fig2 = '_epsilon_crit_vs_gs_Vth';

% out_mat = [outfile, suffix_data, '.mat'];
% out_fig_png = [outfile, suffix_fig, '.png'];
% out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data_gs_pos = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_gs_pos, '.mat']);
outpath_data_pvpw_pos = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_pvpw_pos, '.mat']);
outpath_data_pwinfpv_pos = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_pwinfpv_pos, '.mat']);

outpath_data_gs_original = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_gs_original, '.mat']);
outpath_data_pvpw_original = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_pvpw_original, '.mat']);
outpath_data_pwinfpv_original = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_pwinfpv_original, '.mat']);

outpath_data_epsilon_crit = fullfile(dir_base, dir_data, ...
    [outfile, suffix_data_epsilon, '.mat']);


outpath_fig1_png = fullfile(dir_base, dir_fig, ...
    [outfile, suffix_fig1, '.png']);
outpath_fig1_eps = fullfile(dir_base, dir_fig, ...
    [outfile, suffix_fig1, '.eps']);
% outpath_fig2_png = fullfile(dir_base, dir_fig, ...
%     [outfile, suffix_fig2, '.png']);
% outpath_fig2_eps = fullfile(dir_base, dir_fig, ...
%     [outfile, suffix_fig2, '.eps']);

% data
% save(outpath_data_gs_pos, 'gs_crit_pos')
% save(outpath_data_pvpw_pos, 'pv_pw_pos')
% save(outpath_data_pwinfpv_pos, 'pwinf_pv_pos')
% save(outpath_data_gs_original, 'gs_vth')
% save(outpath_data_pvpw_original, 'pv_pw')
% save(outpath_data_pwinfpv_original, 'pwinf_pv')
% save(outpath_data_epsilon_crit, 'epsilon_crit')
% 
% % figures
% saveas(f1, outpath_fig1_eps,'eps')
% saveas(f1, outpath_fig1_png,'png')
% saveas(f2, outpath_fig2_eps,'eps')
% saveas(f2, outpath_fig2_png,'png')
