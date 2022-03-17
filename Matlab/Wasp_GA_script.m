% Wasp_GA_script.m
% This script runs a Genetic Algorithm (GA) to estimate parameter values 
%  for European wasp (Vespula germanica) using data from Atlas of Living 
%  Australia and a population dynamics model.
% ======================================================================== 
% 
%  The input data is specific to a spatial cluster. In this case we have 
%   two clusters, one for Victoria (VIC) and one for New South Wales (NSW).
%  The results are saved to a vector of GA structures for a number of 
%   iterations using different random seeds.
%               * Written by O.J. Cacho (2021)
%  ========================================================================
clear;
path(path, '.\Data\');
% ------------------load cluster data
%load ga_data_vic
load ga_data_nsw
% ----------------------
parpool('local', 'AttachedFiles',{'Lik_fn'})
%
nruns = 500;
% weights for prediction of [cells invaded, area], and random seed
ga_wt = [0.5, 0.5, 0]; 
%%
for i = 1 : 20
    i
    % ----------------------------- Run GA 
    [GA(i)] = Run_Wasp_GA(nruns,obs_mat, garea, prm.w,  hh_dens, D, ccord, ga_wt);
    ga_wt(3) = ga_wt(3) + 1; % change random seed
   save GA_nsw_55 GA;
end
%
