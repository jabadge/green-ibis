function out = configure_paths(varargin)
%CONFIGURE_PATHS  Show, set, or interactively fill in the toolkit's data paths.
%
%   The MATLAB counterpart of `python -m config.config` on the Python side, writing
%   the same shared config file, so a path set here is also found there.
%
%   Usage:
%      configure_paths()                                % show every setting
%      configure_paths('-set')                          % prompt for missing ones
%      configure_paths('-set','-all')                   % re-prompt for all of them
%      configure_paths('model_dir','/data/models', ...) % set one or more
%      configure_paths('-save')                         % write current paths to disk
%      configure_paths('-clear')                        % drop this session's overrides
%      s = configure_paths('-status')                   % struct array, no printing
%
%   Setting a path with a name/value pair applies for this MATLAB session and
%   takes precedence over the environment and the config file. Follow it with
%   configure_paths('-save') to make it permanent.
%
%   Output:
%      out : for '-save', the config file written (char); for '-status', a
%            struct array with .name/.value/.source; otherwise omitted.
%
%   See also DATA_PATH, PATH_SETTINGS, CONFIG_FILE, PROJECT_CONFIG.

S = path_settings();
args = varargin;

doSet = false; doAll = false; doShow = true;

% --- flags ---------------------------------------------------------------
keep = true(1, numel(args));
for i = 1:numel(args)
    if ~(ischar(args{i}) || isstring(args{i})), continue; end
    switch lower(char(args{i}))
        case '-set',    doSet = true;  keep(i) = false;
        case '-all',    doAll = true;  keep(i) = false;
        case '-clear',  path_state('clear'); keep(i) = false;
        case '-status'
            out = local_status(S); return;
        case '-save'
            out = local_save(S); keep(i) = false; doShow = false;
    end
end
args = args(keep);

% --- name/value pairs ----------------------------------------------------
if mod(numel(args), 2) ~= 0
    error('configure_paths:badArgs', ...
        'Path settings must be given as name/value pairs.');
end
for i = 1:2:numel(args)
    name = char(args{i});
    if ~any(strcmp({S.name}, name))
        error('configure_paths:unknown', ...
            'Unknown path setting "%s"; known: %s.', name, strjoin({S.name}, ', '));
    end
    path_state('set', name, clean_path(args{i+1}));
    doShow = doShow && isempty(args);   % setting explicitly: stay quiet
end
if ~isempty(args), doShow = false; end

% --- interactive fill-in -------------------------------------------------
if doSet
    st = local_status(S);
    for k = 1:numel(S)
        needs = strcmp(st(k).source, 'unset') || ...
                (doAll && ~strcmp(st(k).source, 'built-in'));
        if needs
            try
                data_path(S(k).name);         % prompts and saves
            catch
                % cancelled or not answerable — leave it unset
            end
        end
    end
    doShow = true;
end

% --- report --------------------------------------------------------------
if doShow
    st = local_status(S);
    fprintf('\nData paths   (config file: %s)\n', config_file());
    fprintf('%-18s %-22s %s\n', 'setting', 'source', 'value');
    fprintf('%-18s %-22s %s\n', repmat('-',1,18), repmat('-',1,22), repmat('-',1,40));
    missing = 0;
    for k = 1:numel(st)
        v = st(k).value;
        if isempty(v), v = '(not set)'; missing = missing + 1; end
        fprintf('%-18s %-22s %s\n', st(k).name, st(k).source, v);
    end
    if missing > 0
        fprintf(['\n%d setting(s) not configured — run configure_paths(''-set''), ' ...
                 'or set the GREENIBIS_* environment variables.\n'], missing);
    end
end

end

% ===========================================================================
function st = local_status(S)
%LOCAL_STATUS  Where each setting currently comes from, without prompting.
ov = path_state('get');
[~, data] = config_file();
st = struct('name',{},'value',{},'source',{});
for k = 1:numel(S)
    name = S(k).name;
    value = ''; source = 'unset';
    if isfield(ov, name) && ~isempty(ov.(name))
        value = ov.(name); source = 'session';
    elseif ~isempty(getenv(S(k).env))
        value = clean_path(getenv(S(k).env)); source = S(k).env;
    elseif isfield(data, name) && ~isempty(data.(name))
        value = clean_path(char(data.(name))); source = 'config file';
    else
        try
            value = data_path(name, 'prompt', false); source = 'built-in';
        catch
            value = ''; source = 'unset';
        end
    end
    st(end+1) = struct('name', name, 'value', value, 'source', source); %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------
function f = local_save(S)
%LOCAL_SAVE  Write the currently set paths to the shared config file.
%   Entries for settings that are not set are left as they were, so saving
%   never drops a path another session (or the Python side) put there.
[f, data] = config_file();
st = local_status(S);
for k = 1:numel(st)
    if ~isempty(st(k).value) && any(strcmp(st(k).source, {'session','config file'}))
        data.(st(k).name) = st(k).value;
    end
end
d = fileparts(f);
if ~isfolder(d), mkdir(d); end
fid = fopen(f, 'w');
if fid < 0, error('configure_paths:saveFailed', 'Could not write %s', f); end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(data, 'PrettyPrint', true));
end
