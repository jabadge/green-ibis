function [tnew,rate] = mass2rate(t,mass)
%MASS2RATE: This function converts an ice mass time series to a rate of
%mass change time series
%           rate = \Delta mass / \Delta time
%   INPUTS: 
%           t: vector of times in decimal years
%           mass: vector of ice masses (Gt or some other mass unit)
%   OUTPUTS: 
%           tnew: vector of new times centered on rate time frame
%           rate: mass change rate (e.g., Gt/yr)
%   EXAMPLE: 
%           mass = vol2mass(,...
%                   [md.results.TransientSolution(:).IceVolumeAboveFloatationScaled], ...
%                   md.materials.rho_ice);
%           [tnew, rate] = mass2rate(...
%                               [md.results.TransientSolution(:).time], ...
%                               mass);

rate = diff(mass) ./ diff(t);
tnew = t(1:end-1) + diff(t)/2;
end