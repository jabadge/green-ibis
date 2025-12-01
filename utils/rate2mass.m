function [tnew,mass] = rate2mass(t,rate)
%RATE2MASS: This function converts a time series of rates of mass change 
%into a time series of ice mass change
%
%   INPUTS: 
%           t: vector of times in decimal years
%           rate: vector of rates of ice mass change (e.g., Gt/yr)
%   OUTPUTS: 
%           tnew: vector of new times centered on rate time frame
%           mass: vector of ice mass change (Gt or some other mass unit) 
%   EXAMPLE: 
%           [tnew,mass] = rate2mass(t,rate);

t_temp = diff(t)/2;
tnew = [t(1)-t_temp(1), t(2:end)-t_temp, t(end)+t_temp(end)];
mass = [0 cumsum(rate .* diff(tnew))];
end