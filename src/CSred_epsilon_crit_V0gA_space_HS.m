% Title: Script for detection the minimal epsilon_crit in reduced 
% Connor-Stevens model. This script is specifically for epsilon_crit 
% values above the transcritical (TC) bifn line. 
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

%% Parameters

%Reversal potentials
ENa = 55;
EK = -75;
El = -17;

% conductances
gl = 0.3;
gNa = 120;
gKDR = 20;

% Bifurcation parameters
gA_range = 0:1:120; %medium range
v0_range = -46:1e-1:-25; %medium range
% gA_range = 0:1e-1:120; %large range
% v0_range = -46:1e-2:-25; %large range


options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% Differentials

%with V0 as dynamical variable and without timescale
dVdVC_V0 = @(V, n, gA) -120/((exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^3) - 20*n^4 - gA/((exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^3*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)) - (202661983231672320*exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509)*(V - 55))/(5329489831328509*(exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^4) - (105553116266496*gA*exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013)*(V + 75))/(1820115906168013*(exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^4*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)) - 3/10;
dVdnC_V0 = @(V, n, gA)(73614711755563440*exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371)*(V - 55))/(260839275764371*n^2*(exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)^2*(1/n - 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^3) - 80*n^3*(V + 75) + (4600919484722715*gA*exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963)*(V + 75))/(3215394367055963*n^2*(exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^3*(1/n - 1)*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)^2);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);

%including timescale with V0 as dynamical variable
dVdV_V0 = @(V, n, gA, C) -(120/((exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^3) + 20*n^4 + gA/((exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^3*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)) + (202661983231672320*exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509)*(V - 55))/(5329489831328509*(exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^4) + (105553116266496*gA*exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013)*(V + 75))/(1820115906168013*(exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^4*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)) + 3/10)/C;
dVdn_V0 = @(V, n, gA, C) ((73614711755563440*exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371)*(V - 55))/(260839275764371*n^2*(exp(224164298057848/1304196378821855 - (613455931296362*log(1/n - 1))/260839275764371) + 1)^2*(1/n - 1)*(exp(- (562949953421312*V)/5329489831328509 - 16739365515052668/5329489831328509) + 1)^3) - 80*n^3*(V + 75) + (4600919484722715*gA*exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963)*(V + 75))/(3215394367055963*n^2*(exp(- (35184372088832*V)/1820115906168013 - 2678129957135078/1820115906168013) + 1)^3*(1/n - 1)*(exp(9399834565864460/3215394367055963 - (4600919484722715*log(1/n - 1))/3215394367055963) + 1)^2))/C;
dndV_V0 = @(V, n, V0) (n - 1/(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1))*(exp(- V/80 - 557/800)/320 + 1/(50*(exp(- V/10 - 457/100) - 1)) + (exp(- V/10 - 457/100)*(V/50 + 457/500))/(10*(exp(- V/10 - 457/100) - 1)^2)) + (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715)*(exp(- V/80 - 557/800)/4 - (V/50 + 457/500)/(exp(- V/10 - 457/100) - 1)))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndn_V0 = @(V) (V/50 + 457/500)/(exp(- V/10 - 457/100) - 1) - exp(- V/80 - 557/800)/4;


%% Get TC line
gA_TC = functions.CSred_TC_line(v0_range, options);


%% indices for parfor
v0_idx = 1:length(v0_range);
gA_idx = 1:length(gA_range);

all_gA = ones(length(v0_range),1) * gA_range;
all_v0 = (ones(length(gA_range), 1) * v0_range)';
parfor_paired_values = [all_gA(:) all_v0(:)];


%% Loop

% pre-allocate
epsilon_crit = zeros(length(v0_range), length(gA_range));
gs_crit = zeros(length(v0_range), length(gA_range));
C_crit = zeros(length(v0_range), length(gA_range));

tic;
parfor idx = 1:size(parfor_paired_values,1)
    disp(idx)
    
    % Values for calculations
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    gA = vals(1);
    
    %% (1.) Get SNIC bifurcation info
    [v_sn, n_sn, I_app_sn] = functions.CSred_genV0_sn_coordinates(gA, v0, options);
    
    if ~isnan(n_sn)
        
        %% (2.) Get timescale where left fixed point is Hopf at v_sn
        C_crit(idx) = fsolve(@(C) dVdVC_V0(v_sn, n_sn, gA) ...
            - (C/functions.tau_n(v_sn)), 0, options);
        epsilon_crit(idx) = C_crit(idx)/functions.tau_n(v_sn);
        
        %% (3.) Find gs value for which left fixed point is SNIC at v_sn
        gs_crit(idx) = fsolve(@(gs) ((dVdVC_V0(v_sn, functions.ninf_genV0(v_sn, v0), gA) ...
            * dndnT_V0) + gs), 0, options); %gs is negative, therefore it is '+' instead of '-'
        
        %% (4.) Only collect SN on upper branch and when right fixed point is unstable
        I_app_rfp = I_app_sn + 1e-4;
        
        steady_state = @(y) (-gNa.*functions.minf_gen(y(1)).^3 ...
            .* functions.hinf_gen(functions.inv_ninf_gen(functions.ninf_genV0(y(1), v0))) ...
            * (y(1)-ENa) -gKDR.*functions.ninf_genV0(y(1), v0).^4*(y(1)-EK) ...
            - gA.*functions.mAinf_gen(y(1)).^3 ...
            .* functions.hAinf_gen(functions.inv_ninf_gen(functions.ninf_genV0(y(1), v0))) ...
            * (y(1)-EK) - gl*(y(1)-El) + I_app_rfp);
        
        fixed_pts = functions.CSred_fixed_points_multi(steady_state, options);
        v_rfp = max(fixed_pts);
        n_rfp = functions.ninf_genV0(v_rfp, v0);
        
        % check rfp stability
        jacobian_rfp = [ [dVdV_V0(v_rfp, n_rfp, gA, ...
            C_crit(idx)), dVdn_V0(v_rfp, n_rfp, gA, ...
            C_crit(idx))]; [dndV_V0(v_rfp, n_rfp, v0), ...
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
colormap(jet)
% f1 = imagesc(gA_range, v0_range, log10(epsilon_crit));
f1 = imagesc(gA_range, v0_range, epsilon_crit);
hold on

% colour nan values (rfp exclusion)
set(f1,'alphadata',~isnan(epsilon_crit)) %set nan values to transparent
set(gca,'color','black') %make the background black
set(gcf,'inverthardcopy','off'); 
hold on

% plot transcritical bifn line
plot(gA_TC, v0_range, 'w', 'LineWidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$g_{A}$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} reduced CS model', '\newline', ...
    'Hopf vs SNIC'])


%% Plot gs(V_th) vs epsilon_crit
% f2 = figure;
% for n = 1:size(epsilon_crit, 1)
%     for m = 1:size(epsilon_crit, 2)
%         plot(gs_crit(n,m), epsilon_crit(n,m), 'ok')
%         hold on
% %         plot([0, 0.3], [0, 0.3], 'g-', 'Linewidth', 1.5)
%     end
% end
% % xlim([0, 1.4])
% xlabel('gs'); ylabel('\epsilon');
% title('\epsilon_{crit} vs g_s')



%% Export/Save
outfile = 'CSred_genV0_epsilon_crit_V0gA_space_HS_large_range';

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
% save(outpath_data_eps, 'epsilon_crit')
% save(outpath_data_gs, 'gs_crit')
% save(outpath_data_C, 'C_crit')

% figures
%saveas(f1, out_fig1_eps,'eps')
% saveas(f1, outpath_fig1_png,'png')
%saveas(f2, out_fig2_eps,'eps')
% saveas(f2, outpath_fig2_png,'png')