function cfg = project_config(varargin)
%PROJECT_CONFIG  Central configuration (conventions and known paths) for the toolkit.
%
%   No path here is hard-coded to one machine. The projection/grid conventions
%   are returned directly; data locations are resolved on demand by DATA_PATH
%   (override -> environment variable -> shared config file -> prompt), and the
%   ones already known are copied onto the returned struct for convenience.
%
%   Because a location may legitimately be unconfigured, ALWAYS reach for a
%   data path with DATA_PATH rather than assuming the field exists:
%
%      d = data_path('model_dir');        % resolves, asking once if it must
%      cfg = project_config();  cfg.epsg     % conventions are always present
%
%   Usage:
%      cfg = project_config()
%      cfg = project_config('model_dir','/data/models')   % set, then return
%
%   Name/value options:
%      any setting from PATH_SETTINGS — sets it for this session
%      (equivalent to CONFIGURE_PATHS(name, value))
%
%   Output (struct cfg):
%      .epsg      (double)  projection EPSG code (3413 = Greenland Polar Stereographic)
%      .grid_res  (double)  default regular-grid spacing [m]
%      .obs_bands (cellstr) component order in the Sentinel-type netCDF velocity var
%      plus one field per CONFIGURED path (issm_dir, model_dir, obs_flowline_dir,
%      obs_netcdf_dir, flowline_shp_dir, termini_shp, icefronts_nc, roi_dir) — absent when a
%      path has not been configured yet.
%
%   See also DATA_PATH, CONFIGURE_PATHS, PATH_SETTINGS, SETUP_PATHS.

if ~isempty(varargin)
    configure_paths(varargin{:});
end

% --- conventions ---------------------------------------------------------
cfg.epsg      = 3413;   % NSIDC Sea Ice Polar Stereographic North (Greenland)
cfg.grid_res  = 200;    % [m] default target grid spacing
cfg.obs_bands = {'vv','vx','vy','ex','ey','dT'};  % component dim order in the netCDF

% --- whatever paths are already configured (never prompts) ---------------
S = path_settings();
for k = 1:numel(S)
    try
        cfg.(S(k).name) = data_path(S(k).name, 'prompt', false);
    catch
        % not configured yet: leave the field out, so a typo in a caller
        % surfaces as a clear "no such field" instead of a silent empty path
    end
end

end
