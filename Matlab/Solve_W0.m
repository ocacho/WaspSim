function W0 = Solve_W0(ndetect, delta, hh_dens)
% W0 = Solve_W0(ndetect, delta, hh_dens)
% This function estimates the initial invasion by solving the inverse of  
%   the detection function based on number of detections ndetect. 
% ========================================================================
%   INPUTS:
%     ndetect: number of detections per cell (nc, 1)
%     delta: detection parameter
%     hh_dens: household density per cell (nc, 1)
%   OUTPUT:
%     W0: Matrix of initial wasp invasion (nc, nr) for each cell on the map 
%         and stochastic run.
%   Note: nc is the number of grid cells on the map.
%
%                Written by O.J. Cacho (2021)
% ========================================================================
%
options = optimset('Display', 'off');
u0 = ndetect;
W0 = fsolve(@ObjW,u0,options);

    function xsol = ObjW(x)
      xsol = ndetect - x .* (1 - exp(-delta .* x .* hh_dens)); 
    end

end