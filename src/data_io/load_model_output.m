function M = load_model_output(matpath, varargin)
%LOAD_MODEL_OUTPUT  Load an ISSM model .mat and extract a lean velocity struct.
%
%   Reads a saved ISSM `model` object and pulls out only what is needed for
%   model/observation comparison: the 2-D mesh and the velocity (and optional
%   geometry) time series. Works for both transient runs (TransientSolution)
%   and snapshot inversions (StressbalanceSolution).
%
%   Usage:
%      M = load_model_output(matpath)
%      M = load_model_output(matpath, 'solution','auto', 'trange',[2018 2020])
%
%   Input:
%      matpath : (char) full path OR bare filename found in project_config().model_dir
%   Name/value options:
%      'solution' : (char) 'auto'(default) | 'transient' | 'stressbalance'
%      'trange'   : (1x2 double) [t0 t1] decimal-year window to keep (transient only)
%
%   Output (struct M):
%      .file        (char)              source file
%      .x, .y       (nv x 1 double)     mesh vertex coordinates [m], EPSG per config
%      .elements    (ne x 3 double)     triangulation (1-based vertex indices)
%      .time        (1 x nt double)     decimal years (NaN for a single snapshot)
%      .vx,.vy,.vel (nv x nt double)    velocity components / magnitude [m/yr]
%      .surface     (nv x nt double)    ice surface elevation [m] (if available)
%      .bed         (nv x 1 double)     bed elevation [m] (if available)
%      .solution    (char)              which results field was used

p = inputParser;
p.addParameter('solution', 'auto', @(s) ischar(s) || isstring(s));
p.addParameter('trange', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
p.parse(varargin{:});
opt = p.Results;

% Resolve path (allow bare filename relative to the configured model_dir)
if ~isfile(matpath)
    cfg = project_config();
    cand = fullfile(data_path('model_dir'), matpath);
    if isfile(cand)
        matpath = cand;
    else
        error('load_model_output:notFound', 'Model file not found: %s', matpath);
    end
end

fprintf('Loading model (may be large): %s\n', matpath);
md = loadmodel(matpath);

M = struct();
M.file     = matpath;
M.x        = md.mesh.x(:);
M.y        = md.mesh.y(:);
M.elements = md.mesh.elements;

% Geometry (static) if present
M.bed = [];
try, M.bed = md.geometry.bed(:); catch, end

% Decide which results field to use
res = md.results;
hasT = isfield(res, 'TransientSolution')    && ~isempty(res.TransientSolution);
hasS = isfield(res, 'StressbalanceSolution')&& ~isempty(res.StressbalanceSolution);

sol = lower(char(opt.solution));
if strcmp(sol,'auto')
    if hasT, sol = 'transient'; elseif hasS, sol = 'stressbalance'; else
        error('load_model_output:noResults', ...
            'No TransientSolution or StressbalanceSolution in %s', matpath);
    end
end

switch sol
    case 'transient'
        if ~hasT, error('load_model_output:noTransient','No TransientSolution present.'); end
        S = res.TransientSolution;
        M.time    = [S(:).time];
        M.vx      = [S(:).Vx];
        M.vy      = [S(:).Vy];
        M.vel     = [S(:).Vel];
        M.surface = local_maybe_cat(S, 'Surface');
        % Optional time-window subset
        if ~isempty(opt.trange)
            keep = M.time >= opt.trange(1) & M.time <= opt.trange(2);
            M.time = M.time(keep);
            M.vx = M.vx(:,keep); M.vy = M.vy(:,keep); M.vel = M.vel(:,keep);
            if ~isempty(M.surface), M.surface = M.surface(:,keep); end
        end
        M.solution = 'TransientSolution';

    case 'stressbalance'
        if ~hasS, error('load_model_output:noSnapshot','No StressbalanceSolution present.'); end
        S = res.StressbalanceSolution;
        M.time    = NaN;
        M.vx      = S.Vx(:);
        M.vy      = S.Vy(:);
        M.vel     = S.Vel(:);
        M.surface = local_field_or_empty(md, 'geometry', 'surface');
        M.solution = 'StressbalanceSolution';

    otherwise
        error('load_model_output:badSolution', 'Unknown solution type: %s', sol);
end

fprintf('  mesh: %d vertices, %d elements | %d time step(s) | %s\n', ...
    numel(M.x), size(M.elements,1), numel(M.time), M.solution);

end

% ---------------------------------------------------------------------------
function A = local_maybe_cat(S, fname)
%LOCAL_MAYBE_CAT  Concatenate a per-step field across a solution array, or [].
if isfield(S, fname)
    A = [S(:).(fname)];
else
    A = [];
end
end

function v = local_field_or_empty(md, grp, fname)
%LOCAL_FIELD_OR_EMPTY  Safely fetch md.(grp).(fname) as a column, or [].
v = [];
try
    v = md.(grp).(fname);
    v = v(:);
catch
end
end
