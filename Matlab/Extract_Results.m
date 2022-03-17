function [b_p, w_p, dec_var, prob, ev, sd] = Extract_Results(prm, br_rule, WM, BM, KM, CM)
% [b_p, w_p, dec_var, prob, ev, sd] 
%        = Extract_Results(prm_clust, bioc_rule, WM, BM, KM, CM)
% This function packages the results of a stochastic wasp simulation to 
%   save memory. Results correspond to one experimental treatment.
% ========================================================================
%   INPUTS:
%     prm: parameter struct for wasp and biocontrol *
%     br_rule: biocontrol release rule with fields x_c (coverage), 
%              x_p (intensity per site), x_s (spatial split)  
%     WM: wasp simulation results (nc, nt, nr)    
%     BM: biocontrol simulation results (nc, nt, nr)
%     KM: wasps killed simulation results (nc, nt, nr) 
%     CM: initial biocontrol release (nc, nr)
%     * paramater struct contains:
%       prm.w: wasp paramater struct with fields: alp, kap, mu, gam, h_suit 
%              and delta
%       prm.b: biocontrol paramater struct with fields: alp, kap, mu, gam, 
%              kw and kmu
%   OUTPUTS:
%     b_p: Experimental design vector (1,4) containing biocontrol 
%          parameter set [alp_b, mu_b, gam_b, rho_b]
%     w_p: Parameters for the pest wasp in the experimental treatment
%           [alp_w, kap_w, gam_w, delt_w]    
%     dec_var: Decision variable vector containing the release strategy 
%              (xc, xp) for the treatment in b_p   
%     prob: Probability that a cell will be occupied based on results 
%           from nr stochastic simulations (struct)**  
%     ev: Expected values resulting from nr stochastic simulations 
%         (struct)**  
%     sd: Standard deviations of the results from nr stochastic simulations
%         (struct)** 
%     ** structs contain 4 matrices:
%       W: wasp distribution results (nc, nt)  
%       B: biocontrol distribution results (nc, nt)  
%       K: distribution of wasp killed (nc, nt)  
%       C: distribution of biocontrol cocoon release(nc)  
%
%      The dimensions in the model are:
%       nc: number of grid cells on the map
%       nt: time horizon of simulation
%       nr: number of stochastic (Monte Carlo) runs
%
%                Written by O.J. Cacho (2021)
% ========================================================================
    b_p = zeros(1,4); % biocontrol parameters
    w_p = zeros(1,5); % wasp parameters
    dec_var = zeros(1,3); % decision variables
   % 
    nruns = size(CM,2);
    b_p(1) = prm.b.alp;
    b_p(2) = prm.b.mu ;
    b_p(3) = prm.b.gam;
    b_p(4) = prm.b.kw; % effectiveness (rho_b = phi_b)
    %
    % wasp parameters
    w_p(1) = prm.w.alp;
    w_p(2) = prm.w.kap;
    w_p(3) = prm.w.mu;
    w_p(4) = prm.w.gam;
    w_p(5) = prm.w.delta;
    %
    % decision variables
    dec_var(1) = br_rule.x_c;
    dec_var(2) = br_rule.x_p;
    dec_var(3) = br_rule.x_s;
    %
    % Probability of infestation in cell, yr
    prob.W = sum(WM>0,3)./nruns;
    prob.B = sum(BM>0,3)./nruns;
    prob.K = sum(KM>0,3)./nruns;
    prob.C = sum(CM>0,2)./nruns;
    % means
    ev.W = mean(WM,3);
    ev.K = mean(KM,3);
    ev.B = mean(BM,3);
    ev.C = mean(CM,2);
    % std dev
    sd.W = std(WM,0,3);
    sd.B = std(BM,0,3);
    sd.K = std(KM,0,3);
    sd.C = std(CM,0,2);
end

