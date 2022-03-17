function [n] = Rand_Spread_02(n, K)
% 
% This function Executes one step of random spread for state vector n 
% and probability map K. Each column of K represents the probability that 
% an individual in that column will spread to any row. Columns and rows of 
% square matrix K represent the indexes of n.
% ========================================================================
%   Returns the state vector (n) after random spread based on the 
%   probability map K. The rows of n represent grid cell on the map
%
%   Note: this version of the function perfoms the operation on the whole 
%         matrix at once and can be slow for large maps. In that case use   
%         Rand_Spread.
%
%      Written by O.J. Cacho (2021)
% ========================================================================
%
rspread = max((K-rand(size(K))),0); % random spread matrix
ds = diag(diag(rspread)); % extract diagonal of matrix
rspread = rspread - ds; % keep only nondiagonal elements
ds = eye(size(K)) - diag(sum(rspread)'); % diagonal corrected for spread 
rspread = rspread + ds;
n = rspread * n;
