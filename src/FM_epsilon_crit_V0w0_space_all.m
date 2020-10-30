% Adding together both teps_crit for SNIC/Hopf and teps_crit for SNIC/SN-SH
%
% Author: Jantine Broek
% Created: August 2019
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
v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');


%% Get TC line
w0_TC = functions_FM.FM_TC_line(v0_range, options);


%% Load data
cd(fullfile(dir_base, dir_data));

load('FM_epsilon_crit_V0w0_space_HS_epsilon_crit.mat')
epsilon_crit_HS = teps_crit;
clear teps_crit

load('FM_epsilon_crit_V0w0_space_SNSH_merge.mat');
epsilon_crit_SNSH = epsilon_crit_merge;
clear epsilon_crit_merge

load('FM_epsilon_crit_V0w0_space_HS_rfp.mat')


%% For HS data, set the outer left part to 1e-6 by continuity. 
% The exclusion rfp values are in SNSH data

for i = 1:size(epsilon_crit_HS(:,41:51), 2)
    for j = 1:length(epsilon_crit_HS)
        if epsilon_crit_HS(j,40+i) == 0
            epsilon_crit_HS(j,40+i) = 1e-6;
        end
    end
end

% first three rows are zero or affected by rfp exclusion 
% (will be added again in merge, but exclude for now)
% epsilon_crit_HS(1:3,:) = inf;


%% Add data together

% duplicate matrix
epsilon_crit = epsilon_crit_HS;

for i = 1:size(epsilon_crit,1)
    for j = 1:size(epsilon_crit,2)
        if epsilon_crit(i,j) == 0
            epsilon_crit(i,j) = epsilon_crit_SNSH(i,j);
        elseif isnan(epsilon_crit_SNSH(i,j))
            epsilon_crit(i,j) = NaN;
        end
        
        % select for specific cut-off
        if epsilon_crit(i,j) < 1e-5
            epsilon_crit(i,j) = 1e-5;
        end
    end
end


%% adding exclusion of stable rfp
rfp_coords_add = rfp_coords;

% for now, remove not relevant stable rfp
rfp_coords_add(1:165,:) = 0;
rfp_coords_add(:,1:34) = 0;

% change nans to 0 to add rfp
epsilon_crit(isnan(epsilon_crit)) = 0;

%add rfp exclusion
for i = 1:size(epsilon_crit,1)
    for j = 1:size(epsilon_crit,2)
        
        if rfp_coords_add(i,j) == 1
            epsilon_crit(i,j) = nan;
        else
            epsilon_crit(i,j) = epsilon_crit(i,j);
        end
        
    end
end


%% Change 0 to inf to get a red colour instead of blue

epsilon_crit(epsilon_crit == 0) = Inf;

% remove the blue dots at the lower left corner
epsilon_crit(1:120,92:end) = Inf;
epsilon_crit(16:17,91) = Inf; 
epsilon_crit(27,88) = Inf;
epsilon_crit(30,86) = Inf;

%% Plot data
colormap(jet)
% f1 = imagesc(w0_range, v0_range, log10(epsilon_crit));
f1 = imagesc(w0_range, v0_range, epsilon_crit);

% set rfp stable to black background
set(f1,'alphadata',~isnan(epsilon_crit)) %set nan values to transparent
set(gca,'color','black') %make the background black
set(gcf,'inverthardcopy','off'); 
hold on

%plot transcritical bifn line
plot(w0_TC, v0_range, 'w', 'LineWidth', 2) 

camroll(90)
colorbar

xlabel('$w_0$', 'interpreter','latex','fontsize',12); 
ylabel('$V_0$', 'interpreter','latex','fontsize',12); 
title('\epsilon_{crit} mFHN model')


%% Export/Save

outfile = 'FM_epsilon_crit_V0w0_space_all';

suffix_fig = '';
suffix_data = '';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% figures
% saveas(f1, outpath_fig_eps,'epsc')
% saveas(f1, outpath_fig_png,'png')

