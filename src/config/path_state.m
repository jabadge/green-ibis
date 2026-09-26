function out = path_state(action, name, value)
%PATH_STATE  Session-level path overrides shared by data_path/configure_paths.
%
%   Internal helper: holds the paths set with CONFIGURE_PATHS for the life of the
%   MATLAB session (they take precedence over the environment and the config
%   file, and are not written to disk unless you ask CONFIGURE_PATHS to save).
%
%   Usage:
%      s = path_state('get');
%      s = path_state('set', 'model_dir', '/path/to/models');
%      s = path_state('clear');
%
%   Input:
%      action : (char) 'get' | 'set' | 'clear'
%      name   : (char) setting name (for 'set')
%      value  : (char) path        (for 'set')
%   Output:
%      out : (struct) the current overrides

persistent OVERRIDES
if isempty(OVERRIDES), OVERRIDES = struct(); end

switch lower(action)
    case 'get'
        % nothing to do
    case 'set'
        OVERRIDES.(name) = char(value);
    case 'clear'
        if nargin > 1 && ~isempty(name)
            if isfield(OVERRIDES, name), OVERRIDES = rmfield(OVERRIDES, name); end
        else
            OVERRIDES = struct();
        end
    otherwise
        error('path_state:badAction', 'action must be get, set or clear.');
end
out = OVERRIDES;

end
