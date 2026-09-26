function taub = basal_shear_stress_budd(C,N,ub,p,q)
%BASAL_SHEAR_STRESS_BUDD: This function calculates the basal shear stress
%according to the Budd sliding law. 
%           taub = C^2 * N^(q/p) * ub^(1/p)
%   INPUTS: 
%           C: vector of basal friction coefficients for different
%              locations. Can also be a matrix to include how thickness 
%              evolves in time.
%           N: vector of subglacial effective pressures (Pa). Can also be 
%              a matrix to include how thickness evolves in time.
%           ub: vector of basal sliding velocities (m/s). Can also be 
%              a matrix to include how thickness evolves in time.
%           p: chosen exponent can be value, vector, or matrix
%           q: chosen exponent can be value, vector, or matrix. 
%              Often equal to p.
%   OUTPUTS: 
%           taub: basal shear stress (Pa) vector or matrix depending on 
%                 inputs 
%   EXAMPLE: 
%           taub = basal_shear_stress_budd( ...
%              [md.results.TransientSolution(:).FrictionCoefficient], ...
%              N, ...
%              [md.results.TransientSolution(:).Vel]./md.constants.yts, ...
%              6, ...
%              6);
    
taub = (C.^2).*(N.^(q./p)).*(ub.^(1/p));
end