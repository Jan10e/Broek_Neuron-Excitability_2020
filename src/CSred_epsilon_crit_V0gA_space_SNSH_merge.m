% Title: Merging epsilon critical for SN-SH data in reduced Connor-Stevens
% model. 
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


%% Parameters and functions_CS

% Bifurcation parameters
gA_range = 0:1e-1:120; %large range
v0_range = -46:1e-2:-25; %large range


options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% Load epsilon data
cd(fullfile(dir_base, dir_data));

load('CSred_epsilon_crit_V0gA_space_SNSH_Cstep_1e-3_Cend_1e-2_large_range_epsilon_crit.mat') %smallest range first
epsilon_crit1 = epsilon_crit;
clear epsilon_crit

load('CSred_epsilon_crit_V0gA_space_SNSH_Cstep_1e-2_Cend_1e-1_large_range_epsilon_crit.mat')
epsilon_crit2 = epsilon_crit;
clear epsilon_crit

load('CSred_epsilon_crit_V0gA_space_SNSH_Cstep_1e-1_Cend_3_large_range_epsilon_crit.mat')
epsilon_crit3 = epsilon_crit;
clear epsilon_crit


%% Merge epsilon_crit values
cd(fullfile(dir_base, dir_work));

%merge lowest ranges
epsilon_crit_merge1 = epsilon_crit1; %this should be the smallest range

for i = 1:size(epsilon_crit_merge1,1)
    for j = 1:size(epsilon_crit_merge1,2)
        if isnan(epsilon_crit_merge1(i,j))
            epsilon_crit_merge1(i,j) = epsilon_crit2(i,j);
        else
            epsilon_crit_merge1(i,j) = epsilon_crit1(i,j);
        end
    end
end


%merge all
epsilon_crit_merge2 = epsilon_crit_merge1; %this should be the smallest range

for i = 1:size(epsilon_crit_merge2,1)
    for j = 1:size(epsilon_crit_merge2,2)
        if isnan(epsilon_crit_merge2(i,j))
            epsilon_crit_merge2(i,j) = epsilon_crit3(i,j);
        else
            epsilon_crit_merge2(i,j) = epsilon_crit_merge1(i,j);
        end
    end
end

epsilon_crit_merge = epsilon_crit_merge2;

%% Get TC line
gA_TC = functions_CS.CSred_TC_line(v0_range, options);


%% Plot bifurcation map for
f1 = figure;
colormap(jet)
% imagesc(gA_range, v0_range, log10(epsilon_crit_merge));
imagesc(gA_range, v0_range, epsilon_crit_merge);
hold on

% plot transcritical bifn line
plot(gA_TC, v0_range, 'w', 'LineWidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$g_{A}$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} reduced CS model', '\newline', ...
    'SN-SH'])


%% Export/Save
outfile = 'CSred_epsilon_crit_V0gA_space_SNSH_merge';

suffix_data = '_large_range_epsilon_crit';
suffix_fig = '_large_range';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data, 'epsilon_crit_merge')

% figures
%saveas(f1, outpath_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')



