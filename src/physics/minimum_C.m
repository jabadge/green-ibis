function min_C = minimum_C(CN,H,g,rho_ice)
%MINIMUM_C: This function estimates the minimum basal friction coefficient
%given a CN and ice overburden pressure.
%   INPUTS: 
%           CN: basal friction coefficient * effective pressure. Vector or
%              matrix.
%           H: vector of ice thicknesses (m) for different locations. Can 
%              also be a matrix to include time evolution.
%           g: acceleration due to gravity (m/s) - usually 9.81
%           rho_ice: density of ice (kg/m^3) - usually 917
%   OUTPUTS: 
%           min_C: maximum value of C assuming P_water >= 0
%
%   EXAMPLE: 
%           min_C = minimum_C(CN, ...
%                          [md.results.TransientSolution(:).Thickness], ...
%                          md.constants.g, md.materials.rho_ice);
%   NOTES:
%           N = P_ice - P_water;
min_C = CN ./ (H * g * rho_ice);
end