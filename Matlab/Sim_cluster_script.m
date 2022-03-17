% Sim_cluster_script.m
% This script simulates European Wasp growth in the presence of a 
%  biocontrol agent. Simulations are run for a cluster of grid cells on 
%  a map.
% ========================================================================
% 
%  For details of the model see: 
%   Cacho, O.J. and Hester, S.M. (2022). Modelling biocontrol of invasive 
%         insects: an application to European Wasp (Vespula germanica) in 
%         Australia. Ecological Modeling. 
%         https://doi.org/10.1016/j.ecolmodel.2022.109939
%
%  The model reads .mat files containing maps, model parameters, and 
%   initial conditions and runs a series of stochastic simulations. 
%  The results of each experiment are saved to a different output .mat 
%   file
%  An experiment consists of a set of treatments that represent different 
%   combinations of biocontrol parameters using a full factorial design.
%  Several experiments are run in sequence to evaluate different biocontrol
%   release strategies.
%
%  Simulations are run in parallel for processors with several cores. 
%   A local parallel pool is opened automatically by Matlab when the code 
%   is executed. Alternatively, the user can open the parallel pool 
%   and specify the number of workers using parpool('local', nworkers).
%
%               Written by O.J. Cacho (2021)
%  ========================================================================
clear;
path(path, '.\Data\');
% 
nt = 60; % time periods
%
% (1) =================================== Set biocontrol release strategy 
% one file will be saved for each combination of decision variables
% Coverage: top percent of wasp locations to be inoculated
x_c = [20, 15, 10, 5, 1]; 
% Intensity: proportion of nests to be infested with cocoons per site
x_p = [10, 20, 30, 40, 50]; 
% Split release: 
x_s = [1]; 
%   (1) only top percentage of wasp locations is targeted 
%   (2) top and bottom percentage of wasp locations are targeted (not 
%       functional yet).
%
ff = fullfact([numel(x_c),numel(x_p),numel(x_s)]);
nbr = size(ff,1);
% biocontrol release experiment
br_exp = [x_c(ff(:,1))', x_p(ff(:,2))',x_s(ff(:,3))];
% (2) =========================================== Set experimental design
alpv = [2, 3, 4]; % growth rate
muv = [0.1, 0.3, 0.5]; % mortality 
gamv = [1, 3, 5]; % spread rate
effv = [0.05, 0.1, 0.2]; % effectiveness of biocontrol
%
ff = fullfact([numel(alpv),numel(muv),numel(gamv),numel(effv)]);
nbp = size(ff,1);
% bioc param experiment: 
bp_exp = [alpv(ff(:,1))', muv(ff(:,2))', gamv(ff(:,3))', effv(ff(:,4))'];
%
% ------------------ output file ID based on biocontrol release strategy:
fid = [num2str(br_exp(:,1),'%02d'), repmat('_',nbr,1), ...
       num2str(br_exp(:,2),'%02d')];
%
%% 
% (3) ========================================== Run model for VIC cluster
%
load W0_prm_vic;
fname = [repmat('.\Results\vic_', nbr,1), fid];
% ------------------------------------------ Change effectiveness
for i = 1 : nbr % biocontrol release rules
    % diplay experiment details
    fprintf('\nVIC experiment %d: x_c=%02d, x_p=%02d\n',i, br_exp(i,1), br_exp(i,2));
    % run full experiment
    [b_parm, w_parm, dec_var, prob, ev, sd] = Run_experiment(W0, br_exp(i,:), bp_exp, prm, D, hh_dens, wp_samp, nt);
    save(fname(i,:), 'b_parm', 'w_parm', 'dec_var', 'garea', 'hh_dens', 'prob', 'ev', 'sd');
end
%% 
% (3) ========================================== Run model for NSW cluster
load W0_prm_nsw;
fname = [repmat('.\Results\nsw_', nbr,1), fid];
% ------------------------------------------ Change effectiveness
for i = 1 :  nbr % biocontrol release rules
    % diplay experiment details
    fprintf('\nNSW experiment %d: x_c=%02d, x_p=%02d\n',i, br_exp(i,1), br_exp(i,2));
    % run full experiment
    [b_parm, w_parm, dec_var, prob, ev, sd] = Run_experiment(W0, br_exp(i,:), bp_exp, prm, D, hh_dens, wp_samp, nt);
    save(fname(i,:), 'b_parm', 'w_parm', 'dec_var', 'garea', 'hh_dens', 'prob', 'ev', 'sd');
end
%%
