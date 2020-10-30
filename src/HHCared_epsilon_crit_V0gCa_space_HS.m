% Script for detection the minimal ecrit using the voltage derivative with
% respect to the recovery variable, and the recovery dynamics with respect
% to the voltage
%
% Author: Jantine Broek
% Created: 30 April 2019
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

%% Parameters & Functions

% Reversal potentials
ENa = 55;
EK = -77;
ECa = 85;
El = -54.4;

% Conductances
gl = 0.3;
gNa = 120;
gK = 36;

Ipump = -17;

% Ranges
gCa_range = 0:1e-2:2;
v0_range = -48:5e-2:-29;


options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% Differentials

%with V0 and without timescale
dVdVC_V0 = @(V, n, gCa, V0) (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - gCa*n^3 - 36*n^4 + (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) - 3/10;
dVdnC_V0 = @(V, n, gCa, V0) (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - 144*n^3*(V + 77) - 3*gCa*n^2*(V - 85);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);

%including timescale and with V0 as dynamical variable (to get fixed points)
dVdV_V0 = @(V, n, gCa, C, V0) -(gCa*n^3 - (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 + 36*n^4 - (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) + 3/10)/C;
dVdn_V0 = @(V, n, gCa, C, V0) -(144*n^3*(V + 77) - (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 + 3*gCa*n^2*(V - 85))/C;
dndV_V0 = @(V, n, V0) (n - 1/(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1))*(exp(- V/80 - 3/4)/640 + 1/(100*(exp(- V/10 - 5) - 1)) + (exp(- V/10 - 5)*(V + 50))/(1000*(exp(- V/10 - 5) - 1)^2)) + (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715)*(exp(- V/80 - 3/4)/8 - (V + 50)/(100*(exp(- V/10 - 5) - 1))))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndn_V0 = @(V) (V + 50)/(100*(exp(- V/10 - 5) - 1)) - exp(- V/80 - 3/4)/8;


%% Get TC line
gCa_TC = functions_HH.HHCared_TC_line(v0_range, options);


%% indices for parfor
v0_idx = 1:length(v0_range);
gCa_idx = 1:length(gCa_range);

all_gCa = ones(length(v0_range),1) * gCa_range;
all_v0 = (ones(length(gCa_range), 1) * v0_range)';
parfor_paired_values = [all_gCa(:) all_v0(:)];


%% Loop

% pre-allocate
epsilon_crit = zeros(length(v0_range), length(gCa_range));
gs_crit = zeros(length(v0_range), length(gCa_range));
C_crit = zeros(length(v0_range), length(gCa_range));

tic;
parfor idx = 1:size(parfor_paired_values,1)
    disp(idx)
    
    % Values for calculations
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    gCa = vals(1);
    
    %% (1.) Get SNIC bifurcation info
    [v_sn, n_sn, I_app_sn] = functions.HHCared_sn_coordinates(gCa, v0, options);
    
    if ~isnan(n_sn)
        
        %% (2.) Get timescale where left fixed point is Hopf at v_sn
        C_crit(idx) = fsolve(@(C) dVdVC_V0(v_sn, n_sn, gCa, v0) ...
            - (C/functions.taun60_ab(v_sn)), 0, options);
        epsilon_crit(idx) = C_crit(idx)/functions.taun60_ab(v_sn);
        
        %% (3.) Find gs value for which left fixed point is SNIC at v_sn
        gs_crit(idx) = fsolve(@(gs) ((dVdVC_V0(v_sn, functions.ninf60_genV0(v_sn, v0), ...
            gCa, v0) * dndnT_V0) + gs), 0, options); %gs is negative, therefore it is '+' instead of '-'
        
        %% (4.) Only collect SN on upper branch and when right fixed point is unstable
        I_app_rfp = I_app_sn + 1e-4;
        
        steady_state_rfp =@(y) (-gNa.*functions.minf60_gen(y(1)).^3 ...
            .*(0.89-1.1*functions.ninf60_genV0(y(1), v0)) ...
            *(y(1)-ENa) - gK.*functions.ninf60_genV0(y(1), v0).^4 ...
            *(y(1)-EK) - gl.*(y(1)-El) - gCa.* functions.ninf60_genV0(y(1), v0).^3 ...
            *(y(1)-ECa) + Ipump + I_app_rfp);
        
        fixed_pts = functions.HHCared_fixed_points_multi(steady_state_rfp, options);
        v_rfp = max(fixed_pts);
        n_rfp = functions.ninf60_genV0(v_rfp, v0);
        
        % check rfp stability
        jacobian_rfp = [ [dVdV_V0(v_rfp, n_rfp, gCa, ...
            C_crit(idx), v0), dVdn_V0(v_rfp, n_rfp, gCa, ...
            C_crit(idx), v0)]; [dndV_V0(v_rfp, n_rfp, v0), ...
            dndn_V0(v_rfp)] ];
        
        real_rfp = real(eig(jacobian_rfp));
        
        
        %% (5.) Select for rfp stability
        if epsilon_crit(idx) > 0 && real_rfp(1) > 0
            epsilon_crit(idx) = epsilon_crit(idx);
            gs_crit(idx) = gs_crit(idx);
            C_crit(idx) = C_crit(idx);
        elseif epsilon_crit(idx) > 0 && real_rfp(1) < 0 %rfp stable
            epsilon_crit(idx) = nan;
            gs_crit(idx) = nan;
            C_crit(idx) = nan;
        elseif epsilon_crit(idx) < 0 && real_rfp(1) < 0 %rfp unstable under TC line (rfp is lfp here, due to sign change)
            epsilon_crit(idx) = 0;
            gs_crit(idx) = 0;
            C_crit(idx) = 0;
        elseif epsilon_crit(idx) < 0 && real_rfp(1) > 0 %data under TC line
            epsilon_crit(idx) = nan;
            gs_crit(idx) = nan;
            C_crit(idx) = nan;
        end
        
        
    else %Hopf
        epsilon_crit(idx) = inf;
        gs_crit(idx) = inf;
        C_crit(idx) = inf;
    end
end
toc


%% Bifurcation map
f1 = figure;

subplot(2,1,1)
colormap(jet)
%f1a = imagesc(gCa_range, v0_range, log10(epsilon_crit));
f1a = imagesc(gCa_range, v0_range, epsilon_crit);
hold on

%make NaN black
set(f1a, 'alphadata', ~isnan(epsilon_crit)) %set nan values to transparant
set(gca, 'color', 'black');
set(gcf, 'inverthardcopy', 'off');

% plot transcritical bifn line
plot(gCa_TC, v0_range, 'w', 'LineWidth', 2)

% rotate data
camroll(-90)
set(gca, 'YDir', 'normal')
colorbar
xlabel('$g_{Ca}$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} reduced HH+Ca^+ model', '\newline', ...
    'Hopf vs SNIC'])


subplot(2,1,2)
colormap(jet)
%f1b = imagesc(gCa_range, v0_range, log10(gs_crit));
f1b = imagesc(gCa_range, v0_range, gs_crit);
hold on

% plot transcritical bifn line
plot(gCa_TC, v0_range, 'w', 'LineWidth', 2)

%make NaN black
set(f1b, 'alphadata', ~isnan(gs_crit)) %set nan values to transparant
set(gca, 'color', 'black');
set(gcf, 'inverthardcopy', 'off');

% rotate data
camroll(-90)
set(gca, 'YDir', 'normal')
colorbar
xlabel('$g_{Ca}$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title('$g_s_{crit}', 'interpreter','latex','fontsize',12)



%% Plot gs(V_th) vs epsilon_crit
f2 = figure;
for n = 1:size(epsilon_crit, 1)
    for m = 1:size(epsilon_crit, 2)
        plot(gs_crit(n,m), epsilon_crit(n,m), 'ok')
        hold on
        %         plot([0, 0.3], [0, 0.3], 'g-', 'Linewidth', 1.5)
    end
end
% xlim([0, 1.4])
xlabel('$g_s(V_{th})$', 'interpreter','latex','fontsize',12);
ylabel('$\epsilon_{crit}$', 'interpreter','latex','fontsize',12);
title('\epsilon_{crit} vs g_s(V_{th})')


%% Export/Save
outfile = 'HHCared_epsilon_crit_V0gCa_space_HS';

suffix_data_eps = '_epsilon_crit';
suffix_data_gs = '_gs_crit';
suffix_data_C = '_C_crit';
suffix_fig1 = '';
suffix_fig2 = '_epsilon_crit_vs_gsVth';

out_mat_eps = [outfile, suffix_data_eps, '.mat'];
out_mat_gs = [outfile, suffix_data_gs, '.mat'];
out_mat_C = [outfile, suffix_data_C, '.mat'];
out_fig1_png = [outfile, suffix_fig1, '.png'];
out_fig1_eps = [outfile, suffix_fig1, '.eps'];
out_fig2_png = [outfile, suffix_fig2, '.png'];
out_fig2_eps = [outfile, suffix_fig2, '.eps'];

outpath_data_eps = fullfile(dir_base, dir_data, out_mat_eps);
outpath_data_gs = fullfile(dir_base, dir_data, out_mat_gs);
outpath_data_C = fullfile(dir_base, dir_data, out_mat_C);
outpath_fig1_png = fullfile(dir_base, dir_fig, out_fig1_png);
outpath_fig1_eps = fullfile(dir_base, dir_fig, out_fig1_eps);
outpath_fig2_png = fullfile(dir_base, dir_fig, out_fig2_png);
outpath_fig2_eps = fullfile(dir_base, dir_fig, out_fig2_eps);

% data
save(outpath_data_eps, 'epsilon_crit')
save(outpath_data_gs, 'gs_crit')
save(outpath_data_C, 'C_crit')

% figures
%saveas(f1, out_fig1_eps,'eps')
saveas(f1, outpath_fig1_png,'png')
%saveas(f2, out_fig2_eps,'eps')
saveas(f2, outpath_fig2_png,'png')
