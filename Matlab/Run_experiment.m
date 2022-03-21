function [b_parm, w_parm, dec_var, prob, ev, sd] = Run_experiment(W0, br_exp, bp_exp, prm, D, hh_dens, wp_samp, nt)
% [b_parm, w_parm, dec_var, prob, ev, sd] 
%    = Run_experiment(W0, br_exp, bp_exp, prm, D, hh_dens, wp_samp, nt)
% Runs a simulation experiment for a set of treatments with varying 
%  biocontrol parameters
% ========================================================================
%   INPUTS:
%     W0: Matrix of initial wasp invasion (nc, nr) for each cell on the map 
%         and stochastic run.  
%     br_exp: biocontrol release experiment vector [x_c, x_p, x_s]
%     bp_exp: biocontrol parameter matrix [alp_b, mu_b, gam_b, eps_b] of
%             dimensions (ne,4)
%     prm: parameter struct for wasp and biocontrol * 
%     D: Distance matrix between centroids of each cell and every other 
%        cell on the map 
%     hh_dens: Household density within cells on the map 
%     wp_samp: Matrix of wasp parameters [alp_w, kap_w, gam_w, delt_w] 
%              randomly sampled from GA estimates for nr runs (nr, 54) 
%     nt: number of time periods to simulate (scalar)
%     * parameter struct contains:
%       prm.w: wasp parameter struct with fields: alp, kap, mu, gam, h_suit 
%              and delta
%       prm.b: biocontrol parameter struct with fields: alp, kap, mu, gam, 
%              kw and kmu
%   OUTPUTS:
%     b_parm: Experimental design matrix (ne,4) containing biocontrol 
%          parameter sets [alp_b, mu_b, gam_b, rho_b]
%     w_parm: Wasp parameter matrix [alp_2, kap_w, gam_w, delt_w] 
%             corresponding to rows of b_parm    
%     dec_var: Decision variable matrix containing the release strategy 
%              (xc, xp) for each treatment in b_parm   
%     prob: Probability that a cell will be occupied based on results 
%           from nr stochastic simulations (struct)**  
%     ev: Expected values resulting from nr stochastic simulations 
%         (struct)**  
%     sd: Standard deviations of the results from nr stochastic simulations
%         (struct)** 
%     ** structs contain 4 matrices:
%       W: wasp distribution results (nc, nt, ne)  
%       B: biocontrol distribution results (nc, nt, ne)  
%       K: distribution of wasp killed (nc, nt, ne)  
%       C: distribution of biocontrol cocoon release(nc, ne)  
%
%      The dimensions of matrices in the model are:
%       nc: number of grid cells on the map
%       nt: time horizon of simulation
%       ne: number of experimental treatments
%
%                Written by O.J. Cacho (2021)
% ========================================================================
%
nbp = size(bp_exp,1); % number of biocontrol parameter treatments
b_parm = zeros(nbp,4); % biocontrol [alp, mu, gam, kw]
w_parm = zeros(nbp,5); % wasp [alp, kap, mu, gam, delta]
dec_var = zeros(nbp,3); % decision vars [strat, pct, prop];
% Set biocontrol release rules
bioc_rule.x_c = br_exp(1); 
bioc_rule.x_p = br_exp(2); 
bioc_rule.x_s = br_exp(3); 

elapsed = 0;
for j = 1 : nbp % biocontrol parameters
    %--------------- Set parameters for experiment
    prm.b.alp = bp_exp(j,1);
    prm.b.mu = bp_exp(j,2);
    prm.b.gam = bp_exp(j,3);
    prm.b.rho = bp_exp(j,4); % rho_b
    prm.b.phi = bp_exp(j,4); % phi_b
    %
    tic
    [WM KM, BM CM] = Sim_Biocontrol(W0, bioc_rule, prm, D, nt, hh_dens, wp_samp);
    elapsed = elapsed + toc;
    fprintf(' - treatment: %d, elapsed: %3.2f sec\n',j, elapsed);
    %
    [b_p, w_p, d_var, pr1, ev1, sd1] = Extract_Results(prm, bioc_rule, WM, BM, KM, CM);
    %
    b_parm(j,:) = b_p;
    w_parm(j,:) = w_p;
    dec_var(j,:) = d_var;
    %
    prob.W(:,:,j) = pr1.W;
    prob.B(:,:,j) = pr1.B;
    prob.K(:,:,j) = pr1.W;
    prob.C(:,j) = pr1.C;
    %
    ev.W(:,:,j) = ev1.W;
    ev.B(:,:,j) = ev1.B;
    ev.K(:,:,j) = ev1.W;
    ev.C(:,j) = ev1.C;
    %
    sd.W(:,:,j) = sd1.W;
    sd.B(:,:,j) = sd1.B;
    sd.K(:,:,j) = sd1.W;
    sd.C(:,j) = sd1.C;
end
