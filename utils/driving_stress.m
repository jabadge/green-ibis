function taud = driving_stress(H,g,rho_ice,s)
%DRIVING_STRESS: This function calculates the driving stress
%according to: taud = rho_ice * g * H * sin(alpha) 
%                        
%   INPUTS: 
%           H: vector of thickness (m) at different locations. Can also be 
%               a matrix to include time evolution.
%           g: acceleration due to gravity (m/s) - usually 9.81
%           rho_ice: density of ice (kg/m^3) - usually 917
%           s: slope (degrees) of the ice surface
%   OUTPUTS: 
%           taud: driving stress (Pa) vector or matrix depending on 
%                 inputs 
%   EXAMPLE: 
%           [~,~,s] = surace_slope(md,...
%                             [md.results.TransientSolution(:).Surface],...
%                             2);
%           taud = driving_stress( ...
%              [md.results.TransientSolution(:).Thickness], ...
%              md.constants.g, md.materials.rho_ice, s);
taud = rho_ice * g * H .* sind(s);
end