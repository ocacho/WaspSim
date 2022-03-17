function [Ns, Na] = Grow_Ricker(Ns, alp, bet, mu)
%  [Ns, Na] = Grow_Ricker(Ns, alp, bet, mu)
% This function takes one time step of a Ricker model for growth of a
%   population, where the growth cycle starts in spring and ends in winter.
%   winter mortality occurs on the population present in autumn 
% ========================================================================
% INPUTS
%   Ns: Number of individuals in spring, vector for n sites
%   alp: Ricker growth parameter
%   bet = Ricker growth exponent
%   mu = winter mortality
% Outputs
%   Na = number of indivisuals in autumn at t 
%   Ns = number of individuals in spring at t + 1
%
%      Written by O.J. Cacho (2020)
% ========================================================================
%
Na = alp .* Ns .* exp(-bet .* Ns);
Ns = Na .* (1- mu);
%%