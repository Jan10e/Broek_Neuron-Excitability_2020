% Title: Script for detection the minimal epsilon_crit in reduced 
% Connor-Stevens model. 
%
% Author: Jantine Broek
% Created: 30 April 2019
%% Bifurcation map loop
clear
clc
dbstop if error

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020';

dir_work = '/src';
dir_data = '/data';
dir_fig = '/figures';

cd(fullfile(dir_base, dir_work));

%% Parameter & functions_CS

%ranges
gA_range = 0:1e-1:120; %large range
V0_range = -46:1e-2:-25; %large range

options = optimoptions('fsolve','MaxFunEvals',1e6,...
    'MaxIter',1e6,'TolFun',1e-14,'TolX',1e-14);

%% Load data
cd(fullfile(dir_base, dir_data));

load('CSred_epsilon_crit_V0gA_space_HS_large_range_epsilon_crit.mat');
epsilon_crit_HS = epsilon_crit;
clear epsilon_crit

load('CSred_epsilon_crit_V0gA_space_SNSH_merge_large_range_epsilon_crit.mat');
epsilon_crit_SN = epsilon_crit_merge;
clear epsilon_crit_merge


%% Add data together
cd(fullfile(dir_base, dir_work));

% duplicate matrix
epsilon_crit = epsilon_crit_HS;

for i = 1:size(epsilon_crit,1)
    for j = 1:size(epsilon_crit,2) 
        if epsilon_crit(i,j) == 0
            epsilon_crit(i,j) = epsilon_crit_SN(i,j);
        end
    end
end


%% Get TC line
gA_TC = functions_CS.CSred_TC_line(V0_range, options);


%% In SNSH set NaN to Inf
epsilon_crit(isnan(epsilon_crit)) = Inf;


%% Plot data
colormap(jet)
f1 = imagesc(gA_range, V0_range, log10(epsilon_crit));
% f1 = imagesc(gA_range, v0_range, epsilon_crit);
hold on

% colour nan values (rfp exclusion)
set(f1,'alphadata',~isnan(epsilon_crit_HS)) %set nan values to transparent
set(gca,'color','black') %make the background black
set(gcf,'inverthardcopy','off'); 
hold on

% plot transcritical bifn line
plot(gA_TC, V0_range, 'w', 'LineWidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$V_0$', 'interpreter','latex','fontsize',14); 
ylabel('$\bar{g}_A$', 'interpreter','latex','fontsize',14);
title('Bifn map reduced CS model')


%% Export/Save
outfile = 'CSred_genV0_epsilon_crit_V0gA_space_all';

suffix_fig = '_large_range_log';
suffix_data = '_large_range_epsilon_crit';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data, 'epsilon_crit')

% figures
%saveas(f1, outpath_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')

               