% Title: Merging epsilon critical for SN-SH data
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

%% Load epsilon data
cd(fullfile(dir_base, dir_data));

load('FM_epsilon_crit_V0w0_space_SNSH_epsilon_step_1e-3_epsilon_end_4e-1_epsilon_crit.mat')
epsilon_crit_3 = epsilon_crit;
clear epsilon_crit

load('FM_epsilon_crit_V0w0_space_SNSH_epsilon_step_1e-4_epsilon_end_1e-3_epsilon_crit.mat')
epsilon_crit_4 = epsilon_crit;
clear epsilon_crit

load('FM_epsilon_crit_V0w0_space_SNSH_epsilon_step_1e-5_epsilon_end_1e-4_epsilon_crit.mat');
epsilon_crit_5 = epsilon_crit;
clear epsilon_crit

load('FM_epsilon_crit_V0w0_space_HS_rfp.mat')


%% Ranges and TC line for plot
cd(fullfile(dir_base, dir_work));

v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');

w0_TC = functions_FM.FM_TC_line(v0_range, options);


%% Merge epsilon_crit values

epsilon_crit_merge = epsilon_crit_3;

for i = 1:size(epsilon_crit_4,1)
    for j = 1:size(epsilon_crit_4,2)
        if isinf(epsilon_crit_4(i,j)) || isnan(epsilon_crit_4(i,j))
            epsilon_crit_merge(i,j) = epsilon_crit_merge(i,j);
        else
            epsilon_crit_merge(i,j) = epsilon_crit_4(i,j);
        end
    end
end


for i = 1:size(epsilon_crit_5,1)
    for j = 1:size(epsilon_crit_5,2)
        if isinf(epsilon_crit_5(i,j)) || isnan(epsilon_crit_5(i,j))
            epsilon_crit_merge(i,j) = epsilon_crit_merge(i,j);
        else
            epsilon_crit_merge(i,j) = epsilon_crit_5(i,j);
        end
    end
end


%values change to inf due to rfp incorrect
% epsilon_crit_merge(162:171,93:101) = inf;


%% NaN of epsilon_crit-2 to Inf

for i =1:size(epsilon_crit_merge,1)
    for j = 1:size(epsilon_crit_merge,2)
        if isnan(epsilon_crit_merge(i,j))
            epsilon_crit1(i,j) = inf;
        end
    end
end


%% adding exclusion of stable rfp
rfp_coords_add = rfp_coords;

%for now, remove not relevant stable rfp
rfp_coords_add(1:165,:) = 0;
rfp_coords_add(:,1:34) = 0;

%add rfp exclusion
for i = 1:size(epsilon_crit_merge,1)
    for j = 1:size(epsilon_crit_merge,2)
        if rfp_coords_add(i,j) == 1
            epsilon_crit_merge(i,j) = nan;
        else
            epsilon_crit_merge(i,j) = epsilon_crit_merge(i,j);
        end
    end
end


%% Plot bifurcation map for
f1 = figure;
colormap(jet)
imagesc(w0_range, v0_range, log10(epsilon_crit_merge))
% imagesc(w0_range, v0_range, epsilon_crit_merge);
hold on

%plot transcritical bifn line
plot(w0_TC, v0_range, 'w', 'LineWidth', 2.0) 

camroll(90)
colorbar

xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} mFHN model', '\newline', ...
    'SN-SH'])


%% Export/Save
outfile = 'FM_epsilon_crit_V0w0_space_SNSH_merge';

suffix_fig = '';
suffix_data = '';

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
