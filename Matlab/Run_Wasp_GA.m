function [GA] = Run_Wasp_GA(nruns, obs_mat, garea, z,  hh_dens, D, ccord, wt)
% This function runs t=a genetic algorithm (GA) to minimise the likelihood 
%   function Lik_fn for the data given in obs_mat. 
% ========================================================================
% INPUTS
%   nruns: number of stochastic runs per iteration of GA
%   obs_mat: matrix of observed (reported) presence of wasps (nc, nyr)
%   garea: area of grid cells on the map (nc, 1)
%   z : wasp parameters from prm.w to be modified by GA
%   hh_dens: density of household per grid cell (nc, 1)
%   D : distance matrix (nc, nc)
%   ccord: coordinates of centroids map cells (nc, 2)
%   wt: vector of weights to be applied to predictions based on number 
%        of invaded sites, wtr(1), and area invaded, wtr(2); 
%        wtr(3) is used to pass the random seed  
% OUTPUT
%   GA: struct with fields:
%     .pop: final population of ga solutions (50, 4) with columns
%            representing wasp parameters [alp, kap, gam, delta]             
%     .scores: fitness scores associated with rows of pop (50, 1) 
%     .xstar: paramater set with best (lowest) score (1,4) 
%     .fval: likelihood function value assocaited with xstar 
%     .exitflag=1 when average cumulative change in value of the fitness 
%           function over MaxStallGenerations generations is less than 
%           FunctionTolerance, and the constraint violation is less than 
%           ConstraintTolerance. For other result codes see ga help
%     .output: contains output information returned by the ga function, 
%              see ga help for details
%     .elapsed: time taken to solve the ga (seconds)
%     .D: distance matrix used in the simulations undelying the ga
%     .ccord = coordinates of centroids associated with D
%     .lb: lower bound of ga paramaters (1, 4)
%     .ub: upper bound of ga paramaters (1, 4)
%     .options: struct containing optimisation options for ga (see ga help
%              in Matlab
%  
%                Written by O.J. Cacho (2021)
% ========================================================================
if nargin < 7
    wt = [0.5, 0.5];
end
% GA options
options = optimoptions(@ga);
options = optimoptions(options,'PlotFcn',{@gaplotbestf},'Display','iter');
options = optimoptions(options,'UseParallel', true, 'UseVectorized', false);
% Bounds [alp, kap, gam, delta]
lb = [2, 0.4, 1, 0.1]; 
ub = [5, 2.0, 8, 1.0]; 
% Solve GA
tic
[x_star,fval,exitflag,output,population,scores] = ga({@Lik_fn, obs_mat, garea, z, D, nruns, hh_dens, ccord, wt},4,[],[],[],[],lb,ub,[],options);
elapsed = toc
%
% Package GA results
GA.pop = population;
GA.scores = scores;
GA.xstar = x_star;
GA.fval = fval;
GA.exitflag = exitflag;
GA.output = output;
GA.elapsed = elapsed;
GA.D = D;
GA.ccord = ccord;
GA.lb = lb;
GA.ub = ub;
GA.options = options;


