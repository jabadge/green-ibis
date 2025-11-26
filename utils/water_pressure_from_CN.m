function P_water = water_pressure_from_CN(CN,P_ice,C_est)
%WATER_PRESSURE_FROM_CN: This function estimates water pressure from CN
%given an estimate of C.
%   INPUTS: 
%           CN: basal friction coefficient * effective pressure. Vector or
%              matrix.
%           P_ice: ice overburden pressure (Pa). Vector or matrix.
%           C_est: estimate of the basal friction coefficient. Vector or
%              matrix.
%   OUTPUTS: 
%           P_water: water pressure estimate (Pa)
%
%   EXAMPLE: 
%           P_water = water_pressure_from_CN(CN, ...
%                          [md.results.TransientSolution(:).Thickness], ...
%                          md.constants.g, md.materials.rho_ice, ...
%                          C_est);
%   NOTES:
%           N = P_ice - P_water;
P_water = P_ice - (CN ./ C_est);
if any(P_water<0)
    disp('WARNING: this estimate of C gives some negative water pressures');
end
end