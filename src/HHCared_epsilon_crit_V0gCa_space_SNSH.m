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


%% Parameters and Functions

% Reversal potentials
ENa = 55;
EK = -77;
ECa = 85;
El = -54.4;

% Conductances
gl = 0.3;
gNa = 120;
gK = 36;

I_pump = -17;

% Ranges
gCa_range = 0:1e-2:2;
v0_range = -48:5e-2:-29;
% C_range = 1:1e-2:3;

% C_range in steps
% start with step size (not 0), as 0 is wrongly classified as SNIC (no
% trajectory = no spikes = SNIC), and the algorithm will stop
C_step = 1e-3; %C_step = 1e-1 > dt = 1e-3; C_step = 1e-2 > dt = 1e-4; C_step = 1e-3 > dt = 1e-5
C_start = C_step;
% C_end = C_step / 1e-1; % build until 1.0
C_end = 3;
C_range = C_start:C_step:C_end;

options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);


%% Differentials

% with V0 and without timescale
dVdVC_V0 = @(V, n, gCa) (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - gCa*n^3 - 36*n^4 + (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) - 3/10;
dVdnC_V0 = @(V, n, gCa) (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 - 144*n^3*(V + 77) - 3*gCa*n^2*(V - 85);
dndVT_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndnT_V0 = -1;

dninfdV_V0 = @(V, V0) (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);

%including timescale and with V0 as dynamical variable (to get fixed points)
dVdV_V0 = @(V, n, gCa, C) -(gCa*n^3 - (120*((11*n)/10 - 89/100))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 + 36*n^4 - (101330991615836160*exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269)*((11*n)/10 - 89/100)*(V - 55))/(2666087736960269*(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^4) + 3/10)/C;
dVdn_V0 = @(V, n, gCa, C) -(144*n^3*(V + 77) - (132*(V - 55))/(exp(- (281474976710656*V)/2666087736960269 - 9858542094557778/2666087736960269) + 1)^3 + 3*gCa*n^2*(V - 85))/C;
dndV_V0 = @(V, n, V0) (n - 1/(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1))*(exp(- V/80 - 3/4)/640 + 1/(100*(exp(- V/10 - 5) - 1)) + (exp(- V/10 - 5)*(V + 50))/(1000*(exp(- V/10 - 5) - 1)^2)) + (281474976710656*exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715)*(exp(- V/80 - 3/4)/8 - (V + 50)/(100*(exp(- V/10 - 5) - 1))))/(4600919484722715*(exp((281474976710656*V0)/4600919484722715 - (281474976710656*V)/4600919484722715) + 1)^2);
dndn_V0 = @(V) (V + 50)/(100*(exp(- V/10 - 5) - 1)) - exp(- V/80 - 3/4)/8;


%% Euler solve info
% tend = (5/C_step)*2;
tend = 1e2; %this is already enough for C_step 1e-3
dt = C_step * 1e-2;
t = 0:dt:tend;

% initial values
v_init = 0; %start past fixed points
% n_init = functions.ninf60_genV0(v, v0);
n_init = 0.2;


%% Get TC line
gCa_TC = functions.HHCared_TC_line(v0_range, options);


%% indices for parfor
v0_idx = 1:length(v0_range);
gCa_idx = 1:length(gCa_range);

all_gCa = ones(length(v0_range),1) * gCa_range;
all_v0 = (ones(length(gCa_range), 1) * v0_range)';
parfor_paired_values = [all_gCa(:) all_v0(:)];


%% Bifurcation map loop

% pre-allocate
epsilon_crit = zeros(length(v0_range), length(gCa_range));
C_crit = zeros(length(v0_range), length(gCa_range));
C_branch = zeros(length(v0_range), length(gCa_range));
epsilon_branch = zeros(length(v0_range), length(gCa_range));

tic;
parfor idx = 1:size(parfor_paired_values,1)
    disp(idx)
    
    % Values for calculations
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    gCa = vals(1);
    
    %% (1.) Get SNIC bifurcation info
    [v_sn, n_sn, I_app_sn] = functions.HHCared_sn_coordinates(gCa, v0, options);
    
    %w_SN and Iapp_SN
    if ~isnan(n_sn)
        
        %% (2.) Select for SN at lower branch, i.e. teps_crit<0
        C_branch(idx) = fsolve(@(C) dVdVC_V0(v_sn, n_sn, gCa) ...
            - (C/functions.taun60_ab(v_sn)), 0, options);
        epsilon_branch(idx) = C_branch(idx)/functions.taun60_ab(v_sn);
        
        
        %% Look for SN-SH using selection criteria
        if epsilon_branch(idx) < 0
            
            for i = 1:length(C_range)
                C = C_range(i);
                
                % Get Iapp before saddle-node bifurcation
                I_app = I_app_sn - 1e-3;
                [v, n] = functions.HHCared_forward_euler(v0, gCa, ...
                    C, v_sn, I_app, v_init, n_init, t);
                
                % Get peaks
                [v_max, v_idx] = findpeaks(v, 'MinPeakProminence',20, ...
                    'MinPeakHeight', -30, 'MinPeakDistance',500);
                warning('off','signal:findpeaks:largeMinPeakHeight');
                
                if length(v_max) <= 1 %we start with a lot of peaks. The first value for which there are no multiple peaks is where SNIC happens
                    C_crit(idx) = C;
                    epsilon_crit(idx) = C/functions.taun60_ab(v_sn);
                    break
                end
            end
            
        else %here teps_branch> 0 and on upper branch
            C_crit(idx) = inf;
            epsilon_crit(idx) = inf;
        end
        
        
    else %Hopf only
        C_crit(idx) = inf;
        epsilon_crit(idx) = inf;
    end
    
end
% end
toc


%% Plot bifurcation map for
colormap(jet)
f1 = imagesc(flipud(gCa_range), v0_range, log10(epsilon_crit));
% f1 = imagesc(gCa_range, v0_range, epsilon_crit);
hold on

% plot transcritical bifn line
plot(gCa_TC, v0_range, 'w', 'LineWidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$V_0$', 'interpreter','latex','fontsize',12);
ylabel('$g_{Ca}$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} reduced HH+Ca^+ model', '\newline', ...
    'SN-SH'])

%% Export/Save
outfile = ['HHCared_epsilon_crit_V0gCa_space_SNSH_Cstep_', num2str(C_step), ...
    '_Cend_', num2str(C_end)];

suffix_data_eps = '_epsilon_crit';
suffix_data_C = '_C_crit';
suffix_fig = '_log';

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
