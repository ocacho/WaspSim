function [y] = Lik_fn(u, obs_t, garea, z, D, nruns, hh_dens, cent, wtr)
% [y] = Lik_fn(u, obs_t, garea, z, D, nruns, hh_dens, cent, wtr)
% Likelihood function minimised by Run_Wasp_GA to estimate parameter values
%  for the wasp populatin model based on available detection data. 
% ========================================================================
% INPUTS
%  u : wasp parameter vector [alp, kap, gam, delta] 
%  obs_t: matrix of wasp observations per year for each cell on the map
%  garea : area of grid cells 
%  z : wasp parameters from prm.w to be modified by GA
%  D : distance matrix (nc, nc)
%  nruns : number of stochastic runs
%  hh_km2 : vector of households density (nc, 1)
%  cent: coordinates of centroids map cells (nc, 2)
%  wtr = vector of weights to be applied to predictions based on number 
%        of invaded sites, wtr(1), and area invaded, wtr(2); 
%        wtr(3) is used to pass the random seed  
% OUTPUT
%  y : fitness value to be minimised 
%
%                Written by O.J. Cacho (2021)
% ========================================================================
if nargin < 10
    wt = [0.5, 0.5];
    rseed = 0;
else
    wt = wtr(1:2);
    rseed = wtr(3);
end
rng(rseed,'twister'); %
%
% obs_t is the matrix of observed wasps per (cell, year) in ALA data
nt = size(obs_t,2); % time periods to simulate
% centroids of cells in map
xc = cent(:,1);
yc = cent(:,2);
%
iobs = sum(obs_t,2)>0; % wasp has been reported in cell
xobs = xc(iobs);
yobs = yc(iobs);
% get area of the convex hull around all known infested sites
[cobs, aobs] = convhull(xobs,yobs);  
%
z.alp = u(1);
z.kap = u(2);
z.gam = u(3);
z.delta = u(4); % parameter for prob detection
%
ndetect = obs_t(:,1) ./ garea; % detections/ha in year 1
% inferred infestation in yr 1 based on probability of detection and
% household density:
w0 = Solve_W0(ndetect,z.delta,hh_dens); 
% Estimate initial wasp density
apred = zeros(nruns,1); % infested area predicted
npred = apred; % number of predicted sites infested
for i = 1 : nruns % run stochastic loop
    % Wt is matrix of predictions
    [Wt] = Spread_Stoch_control(w0, z, D, nt, hh_dens);
    % find cells predicted to be invaded at any point in time
    ipred = sum(Wt,2)>0; 
    xpred = xc(ipred);
    ypred = yc(ipred);
    % number of cells predicted to be invaded
    npred(i) = sum(ipred); 
    % area of minimum convex hull around predicted infestations
    try
        [cpred, apred(i)] = convhull(xpred,ypred);
    catch
        % warning('Problem using convhull.  Assigning a value of 0.');
        %apred(1) = 0;
    end
%    end
end
% estimate prediction matrix
% number of sites infested - predicted vs observed
dev1 = (sum(iobs) - mean(npred)).^2;
% area infested - predicted vs observed
dev2 = (aobs - mean(apred)).^2;
% Likelihood function using weights. 
y = wt(1) * dev1 + wt(2) * dev2;

