% Title: Calculation of minimal timescale separation factor needed for bifurcation change
% Check stability rfp. The rfp should be unstable to generate desired
% dynamics. If the rfp is stable, should exclude these values.
%
% Author: Jantine Broek
% Created: Sep 2019
%%
clear
clc
dbstop if error

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020';

dir_work = '/src';
dir_data = '/data';
dir_fig = '/figures';

%% functions_FM and Parameters
w_inf1 = @(v) 2./(1 + exp(-5 * (v)));
w_inf3 = @(v, v0, w0) 2./(1 + exp(-5 * (v - v0))) + w0;

v0_range = -1:1e-2:1;
w0_range = 0.5:-1e-2:-0.5;

epsilon_rfp = 1e-5; %this is the max epsilon-value used in the HS-V0w0 loop

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');

%% Get TC coordinates
w0_TC = functions_FM.FM_TC_line(v0_range, options);

%% Indices
v0_idx = 1:length(v0_range);
w0_idx = 1:length(w0_range);

all_w0 = ones(length(v0_range),1) * w0_range;
all_v0 = (ones(length(w0_range), 1) * v0_range)';
parfor_paired_values = [all_w0(:) all_v0(:)];


%% rfp loop
rfp_coords = nan(1,size(parfor_paired_values,1)); %pre-allocate

tic
parfor idx = 1:size(parfor_paired_values,1)
    disp(idx)
    
    % Values for calculations
    vals = parfor_paired_values(idx,:);
    v0 = vals(2);
    w0 = vals(1);
    
    %% (1) Get bifurcation info
    [v_sn, w_sn, I_app_sn] = functions_FM.FM_sn_coordinates(v0, w0, options);
    
    if ~isnan(w_sn) %three fixed points present > SN
        
        %% (2) Get Iapp right after VSN to see whether rfp is stable
        I_app_rfp = I_app_sn + 1e-3;
        
        steady_state = @(y) (y(1) - (1/3) .* y(1).^3 ...
            - ( w_inf1(y(1) - v0) +w0 ).^2 + I_app_rfp);
        
        %% (3) get value of fixed point left to bifurcation
        fixed_pts = functions_FM.FM_fixed_points_multi(steady_state, options);
        v_rfp = max(fixed_pts);
        w_rfp = w_inf3(v_rfp,v0,w0);
        
        %% (4) Stability rfp
        jacobian_rfp=[ [1 - v_rfp^2, -2 * w_rfp]; [(epsilon_rfp*10*exp(5*v0 - 5*v_rfp))/(exp(5*v0 - 5*v_rfp) + 1)^2, -epsilon_rfp] ];
        
        if min(real(eig(jacobian_rfp))) < 0
            rfp_coords(idx) = 1;
        end
        
        
    end
end
toc

%% Plot bifurcation map for
colormap(jet)
f1 = imagesc(w0_range, v0_range, rfp_coords);
hold on
plot(w0_TC, v0_range, 'w') %plot transcritical bifn line
camroll(90)
colorbar
xlabel('$w_0$', 'interpreter','latex','fontsize',12);
ylabel('$V_0$', 'interpreter','latex','fontsize',12);


%% Export/Save
outfile = 'FM_epsilon_crit_V0w0_space_rfp';

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% data
save(outpath_data, 'rfp_coords')

% figures
%saveas(f1, out_fig_eps,'eps')
saveas(f1, outpath_fig_png,'png')
