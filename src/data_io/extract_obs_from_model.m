function [t,V] = extract_obs_from_model(md)
%EXTRACT_OBS_FROM_MODEL: This function extracts time series of observations
%from the model that were used for transient calibration. Right now, only
%does velocity and assumes vx is the first outputdefinition and vy is the
%second outputdefinition.
%           
%   INPUTS: 
%           md: ISSM model object
%   OUTPUTS: 
%           t: times at which the observations occur (vector, decimal
%           years)
%           V: velocity observations (matrix, m/yr)
%   EXAMPLE: 
%           [t,V] = extract_obs_from_model(md);

% may want to make the positions of these observations in the
% outputdfinitions flexible instead of hard coded
vx = md.outputdefinition.definitions{1}.observations;
vy = md.outputdefinition.definitions{2}.observations;

% can also extract masks of no data values from weights, but not ideal. 
% for now, mask out anywhere that equals 0 exactly, which should deal with 
% no-data locations for how I set up my models.
vx(vx==0)=nan;
vy(vy==0)=nan;

% extract time stamps
if sum(vx(end,:)-vy(end,:))==0
    t = vx(end,:);
else
    error('your vx and vy observations have different time stamps');
end

% calculate velocity magnitude
V = md.constants.yts * ((vx(1:end-1,:).^2)+...
                        (vy(1:end-1,:).^2)).^(1/2);
end

