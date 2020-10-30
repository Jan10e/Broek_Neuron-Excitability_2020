% Title: Calculation of minimal timescale separation factor needed to change from
% bistability (saddle-node on homoclinic cycle) to SNIC. This situation
% occurs under the TC bifn line
%
% Author: Jantine Broek
% Created: July 2019
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

% Reversal potentials
ENa = 55;
EK = -75;
El = -17;

% Conductances
gl = 0.3;
gNa = 120;
gKDR = 20;

% Bifurcation parameters
gA_range = 0:1:120;
v0_range = -46:1e-1:-25;
% gA_range = 0:1e-1:120; %large range
% V0_range = -46:1e-2:-25; %large range
C_range = 1e-2:1e-2:3;


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


%% Euler solve info
tend = 5e1; %choose such that the lowest eps-value (i.e. 1e-2), still gives peaks.
dt = 1e-3;
t = 0:dt:tend;

% initial values
v_init = 0;
n_init = 0.2;

%% Get TC line
gA_TC = functions_CS.CSred_TC_line(v0_range, options);


%% Indices for parfor
% v0_idx = 1:length(v0_range);
% gA_idx = 1:length(gA_range);

all_gA = ones(length(v0_range),1) * gA_range;
all_v0 = (ones(length(gA_range), 1) * v0_range)';
parfor_paired_values = [all_gA(:) all_v0(:)];


%% Loop

% pre-allocate
epsilon_crit = zeros(length(v0_range), length(gA_range));
gs_crit = zeros(length(v0_range), length(gA_range));
C_crit = zeros(length(v0_range), length(gA_range));

C_branch = zeros(length(v0_range), length(gA_range));
epsilon_branch = zeros(length(v0_range), length(gA_range));


tic
parfor idx = 1:size(parfor_paired_values,1)
    disp(idx)
    
    % Values for calculations
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    gA = vals(1);
    
    %% (1.) Get bifurcation info
    [v_sn, n_sn, I_app_sn] = functions_CS.CSred_genV0_sn_coordinates(gA, v0, options);
    
    if ~isnan(n_sn)
        
        %% (2.) Select for SN on lower branch
        C_branch(idx) = fsolve(@(C) dVdVC_V0(v_sn, n_sn, gA) - (C/functions_CS.tau_n(v_sn)), 0, options);
        epsilon_branch(idx) = C_branch(idx)/functions_CS.tau_n(v_sn);
        
        %select for teps_branch on lower branch, i.e. teps_branch<0
        if epsilon_branch(idx) < 0
            
            %% (3.) perform simulation to check for peaks
            for r = 1:length(C_range)
                C = C_range(r)
                
                %Get Iapp before saddle-node bifurcation
                I_app = I_app_sn - 1e-3;
                
                [v, n] = functions_CS.CSred_genV0_forward_euler(v0, gA, C, ...
                    v_sn, I_app, v_init, n_init, t);
                
                
                %Get peaks
                [v_max, v_max_idx] = findpeaks(v, 'MinPeakProminence',5);
                warning('off','signal:findpeaks:largeMinPeakHeight');
                
                
                % select for not too much difference between upper and lower branch of V-nc
%                 [n_max, n_idx] = max(n);
%                 n_min = min(n(n_idx:end));
                
%                 if n_min - n(1) < 0.4
                    if length(v_max) <= 1 %we start with a lot of peaks. The first value for which there are no multiple peaks is where SNIC happens
                        C_crit(idx) = C_range(r);
                        epsilon_crit(idx) = C_range(r)/functions_CS.tau_n(v_sn);
                        break
                    end
%                 end
                
            end
            
        else %here teps_branch > 0 and therefore in SNIC/Hopf situation
            C_crit(idx) = inf;
            epsilon_crit(idx) = inf;
        end
        
    else %no VSN point here, i.e. Hopf only situation
        C_crit(idx) = inf;
        epsilon_crit(idx) = inf;
    end
    
end
toc


%% Plot bifurcation map for
colormap(jet)
% f1 = imagesc(gA_range, v0_range, log10(epsilon_crit))
f1 = imagesc(gA_range, v0_range, epsilon_crit);
hold on

% plot transcritical bifn line
plot(w0_TC, v0_range, 'w', 'LineWidth', 2)

camroll(90)
colorbar
xlabel('$V_0$', 'interpreter','latex','fontsize',12);
ylabel('$g_{A}$', 'interpreter','latex','fontsize',12);
title('\epsilon_{crit} reduced CS model', '\newline', ...
    'SN-SH')


%% Export/Save
outfile = 'CSred_genV0_epsilon_crit_V0gA_space_SNSH';

suffix_data_eps = '_epsilon_crit';
suffix_data_C = '_C_crit';
suffix_fig = '';

out_mat_eps = [outfile, suffix_data_eps, '.mat'];
out_mat_C = [outfile, suffix_data_C, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data_eps = fullfile(dir_base, dir_data, out_mat_eps);
outpath_data_C = fullfile(dir_base, dir_data, out_mat_C);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data_eps, 'epsilon_crit')
save(outpath_data_C, 'C_crit')

% figures
%saveas(f1, out_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')
