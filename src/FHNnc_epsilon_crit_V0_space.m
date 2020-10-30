% Title: get epsilon_crit value for different V0 values in FHN with 
% sigmoidal w-nullcline
%
% Author: Jantine Broek
% Created: Oct 2019
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

%% functions_FHNnc & Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));

v0_range = -0.5:1e-2:1;
% v0_range = -0.14;
epsilon = 0:1e-3:0.5;

epsilon_rfp = 1e-3; %this is the max teps-value used in the loop

options = optimoptions('fsolve','Display','none',...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-14,...
    'TolX',1e-14);

%% Bifurcation map loop

epsilon_crit = zeros(1, length(v0_range));

tic
parfor i = 1:length(v0_range)
    v0 = v0_range(i)
    
    %% (1.) Get SN coordinates
    [v_sn, w_sn, I_app_sn] = functions_FHNnc.FHNnc_sn_coordinates(v0, options);
    
    if ~isnan(w_sn) %SN present
        
        % Get Iapp values for situation before SN and rfp
        I_app = I_app_sn - eps;
        I_app_rfp = I_app_sn + 1e-3;
        
        %% (2.) Get right fixed point (rfp)
        steady_state_rfp = @(y) (y(1) - (1/3) .* y(1).^3 ...
            - ( w_inf1(y(1)-v0) ) + I_app_rfp);
        
        fixed_pts_rfp = functions_FHNnc.FHNnc_fixed_points_multi(steady_state_rfp, options);
        v_rfp = max(fixed_pts_rfp);
        
        %stability rfp
        jacobian_rfp = [ [1-v_rfp^2, -1]; ...
            [(10*epsilon_rfp*exp(5*v0 - 5*v_rfp))...
            /(exp(5*v0 - 5*v_rfp) + 1)^2, -epsilon_rfp] ];
        
        if min(real(eig(jacobian_rfp)))>0 && min(imag(eig(jacobian_rfp)))==0 %also select for Im part as rfp will turn into spiral first
            
            %% (3.) Get left fixed point (lfp)
            steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
                - ( w_inf1(y(1)-v0) ) + I_app);
            
            fixed_pts = functions_FHNnc.FHNnc_fixed_points_multi(steady_state, options);
            v_lfp = min(fixed_pts);
            
            %% (4.) Get eigenvalues of lfp to check whether it changed stability
            for j = 1:length(epsilon)
                jacobian = [ [1-v_lfp^2, -1]; ...
                    [(10*epsilon(j)*exp(5*v0 - 5*v_lfp))...
                    /(exp(5*v0 - 5*v_lfp) + 1)^2, -epsilon(j)]];
                
                real_lambda = real(eig(jacobian));
                
                if min(real_lambda)<0 %SNIC
                    epsilon_crit(i) = epsilon(j);
                    break
                end
            end
            
        else %rfp is stable
           epsilon_crit(i) = NaN;
        end
        
    else %Hopf only
        epsilon_crit(i) = inf;
    end
    
end
toc

%% Plot figure
f1 = imagesc(v0_range, 1, epsilon_crit);
% f1 = imagesc(v0_range, 1, log10(epsilon_crit));
colormap(jet);
set(f1,'alphadata',~isnan(epsilon_crit)); % set nan values to transparent
set(gca,'color','black'); % make the background black
set(gcf,'inverthardcopy','off'); 
colorbar;

xlabel('$V_0$', 'interpreter','latex','fontsize',12); 
title('\epsilon_{crit}', 'fontsize', 12);


%% Export/Save

outfile = 'FHNnc_epsilon_crit_V0w0_space';

suffix_data = '_epsilon_crit';
suffix_fig = '';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
% save(outpath_data, 'epsilon_crit')

% figures
% saveas(f1, outpath_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')

