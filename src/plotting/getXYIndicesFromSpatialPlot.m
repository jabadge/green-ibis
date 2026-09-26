function pos = getXYIndicesFromSpatialPlot(x,y)
%GETXYINDICESFROMSPATIALPLOT: This function lets you click on a plot and
%grab the x, y values. It then finds the vector position of the values in
%the inputs that are closest to the chosen map location.
%           
%   INPUTS: 
%           x: vector x locations
%           y: vector y locations
%   OUTPUTS: 
%           pos: position index of closest x,y pair in input vectors
%   EXAMPLE: 
%           pos = getXYIndicesFromSpatialPlot(md.mesh.x,md.mesh.y)

[x_locs,y_locs]=ginput(1);

x_squared_diff = (x - x_locs(1)).^2;
y_squared_diff = (y - y_locs(1)).^2;
xy_dist = sqrt(x_squared_diff + y_squared_diff);
pos = find(xy_dist==min(xy_dist));
end