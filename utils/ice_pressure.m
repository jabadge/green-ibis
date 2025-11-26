function P_ice = ice_pressure(H,g,rho_ice)
%ICE_PRESSURE: This function calculates the ice overburden pressure
%   INPUTS: 
%           H: vector of ice thicknesses (m) for different locations. Can 
%              also be a matrix to include how thickness evolves in time.
%           g: acceleration due to gravity (m/s) - usually 9.81
%           rho_ice: density of ice (kg/m^3) - usually 917
%   OUTPUTS: 
%           P_ice: ice overburden pressure (Pa) with same dimentions as H 
%   EXAMPLE: 
%           P_ice = ice_pressure( ...
%                      [md.results.TransientSolution(:).Thickness], ...
%                      md.constants.g, ...
%                      md.materials.rho_ice);
P_ice = H * g * rho_ice;
end