function taub = basal_shear_stress_alpha(alpha,ub)
%BASAL_SHEAR_STRESS_ALPHA: This function calculates the basal shear stress
%according to any sliding law using the FrictionAlpha2 output from ISSM. 
%                        
%   INPUTS: 
%           alpha: vector of Alpha2 values output from an ISSM model 
%                  simulation (SI). Can also be a matrix to include time
%                  evolution.
%           ub: vector of basal sliding velocities (m/s). Can also be 
%               a matrix to include time evolution.
%   OUTPUTS: 
%           taub: basal shear stress (Pa) vector or matrix depending on 
%                 inputs 
%   EXAMPLE: 
%           taub = basal_shear_stress_alpha( ...
%              [md.results.TransientSolution(:).FrictionAlpha2], ...
%              [md.results.TransientSolution(:).Vel]./md.constants.yts);
taub = alpha .* ub;
end