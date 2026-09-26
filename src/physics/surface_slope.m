function [sx,sy,s] = surface_slope(md,S,smoothing_amount)
%SURFACE_SLOPE: This function calculates the surface slope
%                        
%   INPUTS: 
%           md: ISSM model object
%           S: surface elevation (m) vector or matrix
%           smoothing_amount: integer 0 or greater. 0 is no smoothing with
%           higher numbers given some smoothing.
%   OUTPUTS: 
%           sx: surface slope (degrees) in x-direction
%           sy: surface slope (degrees) in y-direction
%           s: surface slope (degrees) total
%   EXAMPLE: 
%           [sx,sy,s] = surace_slope(md,...
%                            [md.results.TransientSolution(:).Surface],...
%                            2);
s = zeros(size(S));
sx = zeros(size(S));
sy = zeros(size(S));
for ii = 1:size(S,2)
    [sx_elem,sy_elem,s_elem] = slope(md,S(:,ii));
    s(:,ii) = averaging(md,s_elem,smoothing_amount);
    sx(:,ii) = averaging(md,sx_elem,smoothing_amount);
    sy(:,ii) = averaging(md,sy_elem,smoothing_amount);
end
end