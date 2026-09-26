function taub = basal_shear_stress_schoof(C,N,ub,m,Cmax)
%BASAL_SHEAR_STRESS_SCHOOF: This function calculates the basal shear stress
%according to the Schoof sliding law. 
%                        C^2 ub^m                
%      taub =  ____________________________   
%              (1+(C^2/(Cmax N))^1/m ub )^m  
%   INPUTS: 
%           C: vector of basal friction coefficients for different
%              locations. Can also be a matrix to include time evolution.
%           N: vector of subglacial effective pressures (Pa). Can also be 
%              a matrix to include time evolution.
%           ub: vector of basal sliding velocities (m/s). Can also be 
%              a matrix to include how thickness evolves in time.
%           m: chosen exponent can be value, vector, or matrix
%           Cmax: chosen value, vector, or matrix for Iken's bound 
%                 (typically between 0.17 and 0.84) [SI]
%   OUTPUTS: 
%           taub: basal shear stress (Pa) vector or matrix depending on 
%                 inputs 
%   EXAMPLE: 
%           taub = basal_shear_stress_schoof( ...
%              [md.results.TransientSolution(:).FrictionC], ...
%              N, ...
%              [md.results.TransientSolution(:).Vel]./md.constants.yts, ...
%              3, ...
%              0.5);

numerator=(C.^2).*(ub.^(1/m));
CCmaxNm=((C.^2)./(Cmax*N)).^m;
taub = numerator./((1+(CCmaxNm.*ub)).^(1/m));
end