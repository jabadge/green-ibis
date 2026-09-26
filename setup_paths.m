function cfg = setup_paths(varargin)
%SETUP_PATHS  Put ISSM and this repository's src/ toolkit on the MATLAB path.
%
%   Run this once at the start of a session before calling any toolkit function.
%   For the ISSM *binary* distribution the path setup is simply bin/ + lib/
%   (there is no source tree / devpath to source).
%
%   Data locations are NOT hard-coded: each is resolved on first use by
%   DATA_PATH (session override -> environment variable -> shared config file ->
%   interactive prompt). Nothing is asked for here, so a session that only
%   touches observations never has to answer for the model folder. Run
%   CONFIGURE_PATHS('-set') to fill everything in up front.
%
%   ISSM is optional: without it the observation, ROI and flowline tools all
%   work, and only the model-side functions (LOAD_MODEL_OUTPUT and friends) fail.
%
%   Usage:
%      cfg = setup_paths();
%      cfg = setup_paths('quiet',true);    % no summary printout
%
%   Name/value options:
%      'quiet' : (logical) suppress the path summary (default false)
%
%   Output:
%      cfg : the configuration struct from PROJECT_CONFIG
%
%   See also PROJECT_CONFIG, DATA_PATH, CONFIGURE_PATHS.

p = inputParser;
p.addParameter('quiet', false, @(v) islogical(v)||isnumeric(v));
p.parse(varargin{:});
quiet = logical(p.Results.quiet);

% --- this toolkit: src/ (recursively) plus the script folders ---
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));
addpath(fullfile(here, 'examples'), fullfile(here, 'workflows'));
addpath(genpath(fullfile(here, 'tests')));

% --- ISSM binary distribution (optional) ---
issmOk = false;
try
    issm_dir = data_path('issm_dir', 'prompt', false);
    issm_bin = fullfile(issm_dir, 'bin');
    issm_lib = fullfile(issm_dir, 'lib');
    if isfolder(issm_bin)
        addpath(issm_bin);
        if isfolder(issm_lib)
            addpath(issm_lib);   % compiled mex files (InterpFromMeshToGrid_matlab, etc.)
        end
        issmOk = true;
    else
        warning('setup_paths:issmMissing', ...
            ['ISSM bin folder not found at %s — model-side functions will not ' ...
             'work. Fix it with configure_paths(''issm_dir'', ''/path/to/issm'').'], ...
            issm_bin);
    end
catch
    % issm_dir not configured: only worth mentioning in the summary
end

cfg = project_config();

if ~quiet
    fprintf('\nToolkit paths ready.\n');
    if issmOk
        fprintf('  ISSM             = %s\n', issm_dir);
    else
        fprintf('  ISSM             = (not configured — model-side functions disabled)\n');
    end
    S = path_settings();
    unset = {};
    for k = 1:numel(S)
        n = S(k).name;
        if strcmp(n, 'issm_dir'), continue; end
        if isfield(cfg, n)
            fprintf('  %-16s = %s\n', n, cfg.(n));
        else
            unset{end+1} = n; %#ok<AGROW>
        end
    end
    if ~isempty(unset)
        fprintf('  not configured   : %s\n', strjoin(unset, ', '));
        fprintf('  -> run configure_paths(''-set'') to fill these in (asked once, then remembered)\n');
    end
end

end
