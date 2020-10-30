% Title: Calculation of minimal time-scale separation factor needed for 
% bifurcation change from SNIC bifn to Hopf bifn (above TC line).
%
% Author: Jantine Broek
% Created: 10 Dec 2018
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
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));
w_inf3 = @(v, v0, w0) 2./(1 + exp(-5 * (v - v0))) + w0;

v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;
epsilon_range = 0:1e-5:0.4;

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');


%% Get TC coordinates
w0_TC = functions_FM.FM_TC_line(v0_range, options);


%% Bifurcation map loop

tic
for p = 1:length(v0_range)
    v0 = v0_range(p)
    
    for q = 1:length(w0_range)
        w0 = w0_range(q);
        
        %% (1.) Get SN coordinates
        [v_sn, w_sn, I_app_sn] = functions_FM.FM_sn_coordinates(v0, w0, options);
        
        if ~isnan(w_sn) %three fixed points present > SN
            
            %% (2.) Get left fixed points before SN
            I_app = I_app_sn - eps;
            
            steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
                - ( w_inf1(y(1) - v0) +w0 ).^2 + I_app);
            
            fixed_pts = functions_FM.FM_fixed_points_multi(steady_state, options);
            v_lfp = min(fixed_pts);
            w_lfp = w_inf3(v_lfp,v0,w0);
            
            %% (3.) Get real values of lfp
            for r = 1:length(epsilon_range)
                epsilon = epsilon_range(r);
                
                jacobian = [ [1 - v_lfp^2, -2 * w_lfp]; ...
                    [epsilon * (10 * exp(5*v0 - 5*v_lfp) ...
                    / (exp(5*v0 - 5*v_lfp) + 1)^2), -epsilon] ];
                
                real_lambda = real(eig(jacobian));
                
                if min(real_lambda) < 0 %select epsilon values of negative Re(lambda) values: SNIC
                    epsilon_crit(p,q) = epsilon;
                    break
                end
            end
            
            
        else %Hopf only
            epsilon_crit(p,q) = inf;
        end
    end
end
toc

%% Plot bifurcation map for
colormap(jet)
% f1 = imagesc(w0_range, v0_range, log10(epsilon_crit));
f1 = imagesc(w0_range, v0_range, epsilon_crit);
hold on

%plot transcritical bifn line
plot(w0_TC, v0_range, 'w')

camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);


%% Export/Save
outfile = 'FM_epsilon_crit_V0w0_space_HS';

suffix_fig = '';
suffix_data = 'epsilon_crit';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
% save(outpath_data, 'epsilon_crit')

% figures
% saveas(f1, out_fig_eps,'eps')
% saveas(f1, outpath_fig_png,'png')
