% Title: Calculation of \epsilon_crit, which indicates the change in
% excitability type.
%
% Author: Jantine Broek
% Created: July 2019
%%
clear
clc

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020';

dir_work = '/src';
dir_data = '/data';
dir_fig = '/figures';


%% Parameter & Functions

%ranges
gCa_range = 0:1e-2:2;
v0_range = -48:5e-2:-29;

options = optimoptions('fsolve','MaxFunEvals',1e6,...
    'MaxIter',1e6,'TolFun',1e-14,'TolX',1e-14);


%% Load data
cd(fullfile(dir_base, dir_data));

load('HHCared_epsilon_crit_V0gCa_space_HS_epsilon_crit.mat')
epsilon_crit_HS = epsilon_crit;
clear epsilon_crit

load('HHCared_epsilon_crit_V0gCa_space_SNSH_merge_epsilon_crit.mat')
epsilon_crit_SN = epsilon_crit_merge;
clear epsilon_crit_merge


%% Add data together (old)
cd(fullfile(dir_base, dir_work));

% duplicate matrix
epsilon_crit = epsilon_crit_HS;
epsilon_crit(isinf(epsilon_crit)) = 0;

epsilon_crit_SN(epsilon_crit_SN == 0) = inf;

for i = 1:size(epsilon_crit,1)
    for j = 1:size(epsilon_crit,2)
        
        if epsilon_crit(i,j) == 0
            epsilon_crit(i,j) = epsilon_crit_SN(i,j);
        elseif isnan(epsilon_crit_SN(i,j))
            epsilon_crit(i,j) = NaN;
        end
        
        % select for specific cut-off
        if epsilon_crit(i,j) < 19e-2
            epsilon_crit(i,j) = 1e-3;
        end
    end
end


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
gCa_TC = functions_HH.HHCared_TC_line(v0_range, options);


%% Plot data
colormap(jet)
% f1 = imagesc(gCa_range, v0_range, log10(epsilon_crit));   % for log-scale
f1 = imagesc(gCa_range, v0_range, epsilon_crit);            % without log-scale

% colour nan values (rfp exclusion)
set(f1,'alphadata',~isnan(epsilon_crit)) %set nan values to transparent
set(gca,'color','black') %make the background black
set(gcf,'inverthardcopy','off'); 
hold on

% plot transcritical bifn line
plot(gCa_TC, v0_range, 'w', 'Linewidth', 2)

camroll(-90)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$V_0$', 'interpreter','latex','fontsize',14); 
ylabel('$\bar{g}_{Ca}$', 'interpreter','latex','fontsize',14);
title('Bifn map reduced HH+Ca^+ model')


%% Export/Save
outfile = 'HHCared_epsilon_crit_V0gCa_space_all';

suffix_data = '_epsilon_crit';
suffix_fig = '';

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
