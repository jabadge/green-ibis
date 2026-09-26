function CN = CN_from_Budd(taub,ub,p,q)
%CN_FROMS_BUDD: This function calculates CN from the basal shear stress
%according to the Budd sliding law. 
%           (C^2)^(p/q) * N = (taub / ub^(1/p))^(p/q)
%   INPUTS: 
%           taub: vector of basal shear stress (Pa). Can also be 
%              a matrix to include time evolution.          
%           ub: vector of basal sliding velocities (m/s). Can also be 
%              a matrix to include time evolution.
%           p: chosen exponent can be value, vector, or matrix
%           q: chosen exponent can be value, vector, or matrix. 
%              Often equal to p.
%   OUTPUTS: 
%           CN: basal friction coefficient * effective pressure 
%   EXAMPLE: 
%           CN = CN_from_Budd(taub, ...
%              [md.results.TransientSolution(:).Vel]./md.constants.yts, ...
%              6, ...
%              6);
CN = (taub ./ (ub.^(1/p))).^(p/q);
end