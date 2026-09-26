function mass = vol2mass(volume,rho_ice)
%VOLUME2MASS: This function converts ice volume to mass
%           mass = volume * rho_ice / 10^12
%   INPUTS: 
%           volume: scalar, vector, or matrix of volume(s) (m^3)
%           rho_ice: density of ice (kg/m^3) - usually 917
%   OUTPUTS: 
%           mass: ice mass (Gt)
%   EXAMPLES: 
%           mass = vol2mass(,...
%                   [md.results.TransientSolution(:).IceVolumeScaled], ...
%                   md.materials.rho_ice);
%           mass = vol2mass(,...
%                   [md.results.TransientSolution(:).IceVolumeAboveFloatationScaled], ...
%                   md.materials.rho_ice);

mass = volume * rho_ice / 10^12;
end