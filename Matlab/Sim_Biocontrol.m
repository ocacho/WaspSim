function [WM, KM, BM, CM] = Sim_Biocontrol(W0, bioc_rule, prm, D, nt, hh_dens, wp_samp)
% [WM, KM, BM, CM] 
%    = Sim_Biocontrol(W0, bioc_rule, prm, D, nt, hh_dens, wp_samp)
% This function runs a biocontrol simulation using random wasp parameters 
%   from wp_samp.
% ========================================================================
%   INPUTS:
%     W0: Matrix of initial wasp invasion (nc, nr) for each cell on the map 
%         and stochastic run.  
%     bioc_rule: biocontrol release rule with fields x_c (coverage), 
%              x_p (intensity per site), x_s (spatial split)  
%     prm: parameter struct for wasp and biocontrol * 
%     D: Distance matrix between centroids of each cell and every other 
%        cell on the map (nc, nc) 
%     nt: number of time periods to simulate
%     hh_dens: Household density within cells on the map (nc, 1)
%     wp_samp: Matrix of wasp parameters [alp_w, kap_w, gam_w, delt_w] 
%              randomly sampled from GA estimates for nr runs (nr, 4) 
%     * parameter struct contains:
%       prm.w: wasp parameter struct with fields: alp, kap, mu, gam, h_suit 
%              and delta
%       prm.b: biocontrol parameter struct with fields: alp, kap, mu, gam, 
%              kw and kmu
%   OUTPUTS:
%       WM: wasp distribution results (nc, nt, nr)  
%       KM: distribution of wasp killed (nc, nt, nr)  
%       BM: biocontrol distribution results (nc, nt, nr)  
%       CM: distribution of biocontrol cocoon release(nc, nr)  
%    The dimensions of matrices in the model are:
%       nc: number of grid cells on the map
%       nt: time horizon of simulation
%       nr: number of stochastic runs
%
%                Written by O.J. Cacho (2021)
% ========================================================================
%
[ng,nruns] = size(W0); 
% Variables to save results
WM = zeros(ng,nt,nruns); % nests / km2
BM = zeros(ng,nt,nruns); % biocontrol adults / km2
KM = zeros(ng,nt,nruns); % nests killed / km2
CM = zeros(ng,nruns); % initial cocoons
parfor i = 1 : nruns
    rng(i,'twister'); % for replicability
    %
    pw = struct();
    pw = prm.w;
    pw.alp = wp_samp(i,1);
    pw.kap =  wp_samp(i,2);
    pw.gam =  wp_samp(i,3);
    pw.delta =  wp_samp(i,4);

    %
    w0 = W0(:,i);
    c0 = zeros(ng,1); 
    idx = find(w0 >= prctile(w0,(100 - bioc_rule.x_c)));
    c0(idx) = w0(idx) * (bioc_rule.x_p ./ 100); % innoculate x% of nests
    %
    [Wt, Kt, Bt] = Spread_Stoch_Biocontrol(w0, c0, pw, prm.b, D, nt, hh_dens);
    WM(:,:,i) = Wt; 
    KM(:,:,i) = Kt; 
    BM(:,:,i) = Bt; 
    CM(:,i) = c0; 
end
