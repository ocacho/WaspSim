function [Wt, Kt, Bt] = Spread_Stoch_Biocontrol(w0, c0, wp, bp, D, nt, hh_dens)
%
% [Wt, Kt, Bt] = Spread_Stoch_Biocontrol(w0, c0, wp, bp, D, nt, hh_dens)
% This function executes one stochastic run of population growth and spread 
%   of wasps and biocontrol agent for given inital conditions and parameter 
%   values
% ========================================================================
%   INPUTS:
%    w0 : inital number of wasp nests in spring (nc, 1)
%    c0 : inital number of surviving cocoons released (nc, 1)
%    wp : struct of wasp parameters with fields alp, kap, mu, gam, h_suit 
%         and delta
%    bp : struct of biocontrol parameters with fields alp, kap, mu, gam, 
%         kw and kmu
%    D : Distance matrix (nc, nc) 
%    nt : number of years to simulate 
%    hh_dens : households density (nc, 1)
%
%   OUTPUTS:
%    Wt : matrix of wasp nests density (nc, nt)
%    Kt : matrix of nests killed (nc, nt)
%    Bt : matrix of biocontrol density (nc, nt)
%
%    The dimensions of matrices in the model are:
%       nc: number of grid cells on the map
%       nt: time horizon of simulation
%                Written by O.J. Cacho (2021)
% ========================================================================
nw = length(w0);
Wt = zeros(nw,nt); % nests per ha
Kt = zeros(nw,nt); % killed
Bt = zeros(nw,nt); % biocontrol adults per ha
wt = w0;
bt = (1 - bp.mu) .* c0; % survival of cocoons
kv = wp.h_suit .* wp.kap;
bet_v = (log(wp.alp)./kv); % beta implied by habitat suitability
Kw = Disp_Kernel(wp.gam, D); 
Kb = Disp_Kernel(bp.gam, D); 
% 
for t = 1 : nt
    % adjust parameters for wasp-biocontrol interaction
    bet_b = log(bp.alp)./(bp.kap .* wt); % biocontrol beta depends on wasp nests
    bet_b(bet_b==inf) = 1000; % if no wasp nests biocontrol can't establish
    alp_bw = min((bp.kw .* bt),1); % effect of biocontrol on alpha
    alp_w = wp.alp .* (1 - alp_bw); % wasp alpha depends on biocontrol density
    mu_w = min((bp.kmu .* bt), 1); % effect of biocontrol on wasp winter mortality
    %
    wt = Rand_Spread(wt,Kw); % wasp spread occurs first
    bt = Rand_Spread(bt,Kb); % biocontrol spread
    wt = Grow_Ricker(wt, alp_w, bet_v, mu_w); % then population growth on site
    bt = Grow_Ricker(bt, bp.alp, bet_b, bp.mu); % then population growth on site
    % detect and destroy
    p_detect = 1 - exp(- wp.delta .* wt .* hh_dens); % 
    destroy = wt .* p_detect;
    wt = wt - destroy; % all detections are destroyed
    Wt(:,t) = wt;
    Kt(:,t) = destroy;
    Bt(:,t) = bt;
end
