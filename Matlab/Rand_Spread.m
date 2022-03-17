function [n] = Rand_Spread(n, K)
% 
% This function Executes one step of random spread for state vector n 
% and probability map K. Each column of K represents the probability that 
% an individual in that column will spread to any row. Columns and rows of 
% square matrix K represent the indexes of n.
% ========================================================================
%   Returns the state vector (n) after random spread based on the 
%   probability map K. The rows of n represent grid cell on the map.
%
%      Written by O.J. Cacho (2021)
% ========================================================================
%
nr = length(n);
idx = find(n > 0);
ni = length(idx);
n0 = n;
n_delta = zeros(nr,ni); % matrix of changes in numbers
for i = 1 : ni
    j = idx(i);
    k = K(:,j); % kernel column
    n1 = k .* n(j); % spread for column
    ndisp = n0(j) - n1(j);
    is_on = rand(nr,1) < k;
    if sum(is_on) > 0
        n_delta(j) = n_delta(j) - ndisp; % take them away from here
        n_delta(is_on) = n_delta(is_on) + ndisp / sum(is_on); % spread them evenly
    end
end
n = n + sum(n_delta,2);

