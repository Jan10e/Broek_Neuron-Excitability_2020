% Title: Merging epsilon critical for SN-SH data
% This code depends on `HHCared_epsilon_crit_V0gCa_space_SNSH.m'
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

% Ranges
gCa_range = 0:1e-2:2;
v0_range = -48:5e-2:-29;

options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% Load epsilon data
cd(fullfile(dir_base, dir_data));


load('HHCared_epsilon_crit_V0gCa_space_SNSH_Cstep_0.001_Cend_0.01_epsilon_crit.mat') %smallest range first
epsilon_crit1 = epsilon_crit;
clear epsilon_crit

load('HHCared_epsilon_crit_V0gCa_space_SNSH_Cstep_0.01_Cend_0.1_epsilon_crit.mat') 
epsilon_crit2 = epsilon_crit;
clear epsilon_crit

load('HHCared_epsilon_crit_V0gCa_space_SNSH_Cstep_0.1_Cend_3_epsilon_crit.mat')
epsilon_crit3 = epsilon_crit;
clear epsilon_crit


%% Merge epsilon_crit values
cd(fullfile(dir_base, dir_work));

%merge lowest ranges
epsilon_crit_merge1 = epsilon_crit1; %this should be the smallest range

for i = 1:size(epsilon_crit_merge1,1)
    for j = 1:size(epsilon_crit_merge1,2)
        if epsilon_crit_merge1(i,j) == 0
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
        if epsilon_crit_merge2(i,j) == 0
            epsilon_crit_merge2(i,j) = epsilon_crit3(i,j);
        else
            epsilon_crit_merge2(i,j) = epsilon_crit_merge1(i,j);
        end
    end
end

epsilon_crit_merge = epsilon_crit_merge2;



%% Get TC line
gCa_TC = functions_HH.HHCared_TC_line(v0_range, options);


%% Plot bifurcation map for
f1 = figure;
colormap(jet)
% imagesc(gCa_range, v0_range, log10(epsilon_crit_merge));
imagesc(gCa_range, v0_range, epsilon_crit_merge);
hold on

% plot transcritical bifn line
plot(gCa_TC, v0_range, 'w', 'LineWidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$g_{Ca}$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} reduced HH+Ca^+ model', '\newline', ...
    'SN-SH'])


%% Export/Save
outfile = 'HHCared_epsilon_crit_V0gCa_space_SNSH_merge';

suffix_data = '_epsilon_crit';
suffix_fig = '';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data, 'epsilon_crit_merge')

% figures
%saveas(f1, out_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')



