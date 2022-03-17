function [Wt, Kt] = Spread_Stoch_control(w0, z, D, nt, hh_dens)
% [Wt, Kt] = Spread_Stoch_control(w0, z, D, nt, hh_ha)
% This function executes one stochastic run of population growth and spread 
%   of wasps for given inital conditions and parameter values. It assumes
%   that all wasp nests detected are destroyed.
% ========================================================================
%   INPUTS:
%    w0 : inital number of wasp nests in spring (nc, 1)
%    z : struct of wasp parameters with fields alp, kap, mu, gam, h_suit 
%         and delta
%    D : Distance matrix (nc, nc) 
%    nt : number of years to simulate 
%    hh_dens : households density (nc, 1)
%
%   OUTPUTS:
%    Wt : matrix of wasp nest density (nc, nt)
%    Kt : matrix of nests killed (nc, nt)
%
%    The dimensions of matrices in the model are:
%       nc: number of grid cells on the map
%       nt: time horizon of simulation
%
%                Written by O.J. Cacho (2021)
% ========================================================================
%
% variables
nw = length(w0);
Wt = zeros(nw,nt); % nests per ha
Kt = zeros(nw,nt); % killed
wt = w0;
kv = z.h_suit .* z.kap;
bet_v = (log(z.alp)./kv); % beta implied by habitat suitability
K = Disp_Kernel(z.gam, D); 
% 
for t = 1 : nt
    wt = Rand_Spread(wt,K); % spread occurs first
    wt = Grow_Ricker(wt, z.alp, bet_v, z.mu); % then population growth on site
    % detact and destroy
    p_detect = 1 - exp(- z.delta .* wt .* hh_dens); % 
    destroy = wt .* p_detect;
    wt = wt - destroy; % all detections are destroyed
    Wt(:,t) = wt;
    Kt(:,t) = destroy;
end
