function dn = decyear2datenum(yr)
%DECYEAR2DATENUM  Convert decimal year to MATLAB datenum.
%
%   Inverse of DATENUM2DECYEAR. Handy for placing ISSM model times
%   (decimal years) onto a datetime axis alongside observed series.
%
%   Usage:
%      dn = decyear2datenum(yr)
%
%   Input:
%      yr : (numeric array) decimal years
%   Output:
%      dn : (double array, same shape) MATLAB serial date numbers

sz = size(yr);
yr = yr(:);
y  = floor(yr);
frac = yr - y;
yearStart = datenum(y,     1, 1);
yearEnd   = datenum(y + 1, 1, 1);
dn = yearStart + frac .* (yearEnd - yearStart);
dn = reshape(dn, sz);

end
