function C_est = C_estimate_no_water(CN,P_ice)
%C_ESTIMATE_NO_WATER: This function estimates C assuming water pressure
%reaches 0 once during the time series
%   INPUTS: 
%           CN: basal friction coefficient * effective pressure. Vector or
%              matrix.
%           P_ice: ice overburden pressure (Pa). Vector or matrix.
%   OUTPUTS: 
%           C_est: estimate of the basal friction coefficient. Vector or
%              matrix.
%
%   EXAMPLE: 
%           C_est = C_estimate_no_water(CN,P_ice);
C_est = max(CN./P_ice,[],2,'omitnan');
end