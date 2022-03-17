function [K] = Disp_Kernel(gam, D) 
% [K] = Disp_Kernel(gam, D) 
% This function returns a Cauchy dispersal kernel (K) for distance 
%  matrix D and with spread parameter gam the dimensions of K are the same
%  as those of D
% ========================================================================
%  INPUTS:
%  - gam: dispersal parameter (median dispersal distance)  
%  - D: distance matrix or vector of distances between points on the map 
%
%  OUTPUT:
%  - K: probability map normalised so that sum of each column = 1 
%
%               * Written by O.J. Cacho (2020)
% ========================================================================
%
[nr] = size(D,1);
kernel = 1 ./ ((pi .* gam) .* (1 + (D ./ gam).^2)); % Cauchy 
sk = sum(kernel);
K = kernel ./ (repmat(sk,nr,1));% normalise

