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
dir_base = '/Users/jantinebroek/Documents/03_Prog_Comp/01_Theory/03_TypeI_II/FHNm/';

dir_work = '01_code/matlab';
cd(fullfile(dir_base, dir_work));

dir_data = '02_data/data_matlab';
dir_fig = '03_figures/figures_matlab';


%% Functions and Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));
w_inf3 = @(v, v0, w0) 2./(1 + exp(-5 * (v - v0))) + w0;

v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;
epsilon_range = 1e-3:1e-3:4e-1; %tend = 1e4
% epsilon_range = 1e-4:1e-4:1e-3; %tend = 1e5
% epsilon_range = 0.0:1e-5:1e-4; #tend = 1e6

% epsilon_range = 1e-2:1e-2:0.7;
% epsilon_range = 1e-3:1e-3:0.2; %narrow range to look at 1e-3. Changes tend = 1e4
% epsilon_range = 1e-4:1e-4:0.06; %narrow range to 1e-4 under TC line. Change tend = 1e5


% epsilon_rfp = 1e-5; %this is the max teps-value used in the HS-v0w0 loop

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');

%% Euler solve info
tend = 1e4; %choose such that the lowest eps-value (i.e. 1e-2), still gives peaks.
dt = 0.01;
t = 0:dt:tend;

v_init = 1.0;
w_init = 0;

%% Get TC coordinates
w0_TC = functions.FM_TC_line(v0_range, options);

%% Indices
v0_idx = 1:length(v0_range);
w0_idx = 1:length(w0_range);

all_w0 = ones(length(v0_range),1) * w0_range;
all_v0 = (ones(length(w0_range), 1) * v0_range)';
parfor_paired_values = [all_w0(:) all_v0(:)];

%% Bifurcation map loop
epsilon_crit = nan(1,size(parfor_paired_values,1));

parfor idx = 1:size(parfor_paired_values, 1)
    disp(idx)
    
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    w0 = vals(1);
    
    %% (1) Get SN bifurcation info
    [v_sn, w_sn, I_app_sn] = functions.FM_sn_coordinates(v0, w0, options);
    
    %w_SN and I_app_SN
    if ~isnan(w_sn) &&  w_sn <= 0 %w_sn > 0 and therefore in SNIC/Hopf situation
        
        %Get I_app before saddle-node bifurcation
        I_app = I_app_sn - 1e-3;
        %             I_app_rfp = I_app_sn + 1e-3;
        
        for r = 1:length(epsilon_range)
            %             epsilon = epsilon_range(r)
            
            %% (2.) Get rfp stability
            steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
                - ( w_inf1(y(1) - v0) +w0 ).^2 + I_app);
            
            fixed_pts = functions.FM_fixed_points_multi(steady_state, options);
            
            v_rfp = max(fixed_pts);
            w_rfp = w_inf3(v_rfp, v0, w0);
            
            %stability rfp
            jacobian_rfp=[ [1 - v_rfp^2, -2 * w_rfp]; ...
                [(epsilon_range(r)*10*exp(5*v0 - 5*v_rfp))...
                /(exp(5*v0 - 5*v_rfp) + 1)^2, -epsilon_range(r)] ];
            
            %% (3.) Get peaks to determine epsilon_crit for SN-SH to SNIC
            if min(real(eig(jacobian_rfp))) > 0
                
                [v, w] = functions.FM_forward_euler(v0, w0, epsilon_range(r), ...
                    I_app, v_init, w_init, t);
                
                [v_max, ~] = findpeaks(v(1,:), 'MinPeakProminence',0.05, ...
                    'MinPeakHeight',-0.8);
                warning('off','signal:findpeaks:largeMinPeakHeight');
                
                if length(v_max) <= 1 %we start with a lot of peaks. The first value for which there are no multiple peaks is where SNIC happens
                    epsilon_crit(idx) = epsilon_range(r);
                    break
                else
                    epsilon_crit(idx) = NaN; %here it didn't manage to have a eps_crit value before 0.4. Most of the time the rfp will turn stable after 0.4
                end
                
            else %the rfp was stable
                epsilon_crit(idx) = NaN;
            end
        end
        
    else %or Hopf only or w_SN > 0 and therefore in SNIC/Hopf situation
        epsilon_crit(idx) = inf;
    end
    
end



%% Plot bifurcation map for
colormap(jet)
% f1 = imagesc(w0_range, w0_range, log10(epsilon_crit))
% f1 = imagesc(w0_range, v0_range, epsilon_crit);
hold on

%plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar

xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);
title(['\epsilon_{crit} mFHN model', '\newline', ...
    'SN-SH'])


%% Get rid of zeros
% epsilon_critN = epsilon_crit;
% epsilon_critN(find(epsilon_critN == 0)) =inf;


%% Export/Save
outfile = 'FM_epsilon_crit_V0w0_space_SNSH_epsilon_step_1e-4_epsilon_end_1e-3';

suffix_fig = '';
suffix_data = '_epsilon_crit';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data, 'epsilon_crit')

% figures
%saveas(f1, out_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')
