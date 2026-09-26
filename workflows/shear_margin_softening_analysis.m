%SHEAR_MARGIN_SOFTENING_ANALYSIS
%
% Replicates the per-catchment shear margin softening analysis from the
% '<domain>_Rheology' step of runme_greenibis.m, and exports the diagnostic
% fields needed to build a per-catchment summary PDF (see
% plot_shear_margin_summary.py) as plain (non-classdef) .mat files.
%
% Unlike runme_greenibis.m, this script does not modify or re-save a model:
% the softening and resulting stress balance re-solve is already baked into 
% <domain>_Rheology.mat model output. This script only redoes the 
% extraction/correlation logic per catchment (which is not saved anywhere) 
% so that it can label which vertices were softened, and pulls the matching 
% before/after fields out of two saved models:
%
%   model_before    - stress balance solved with the ORIGINAL rheology
%                      (runme_greenibis.m 'Friction' step output)
%   model_after     - same mesh, with rheology_B AFTER shear margin
%                      softening (runme_greenibis.m 'Rheology' step output)
%
% Both must share the same mesh (none of these steps remesh), which
% is true for the standard runme_snapshot.m pipeline.
%
% Usage:
%   Edit the CONFIGURATION block below to match your environment, then:
%       matlab -batch "shear_margin_softening_analysis"
%   or from within MATLAB with ISSM devpath already loaded:
%       shear_margin_softening_analysis
%
% This script depends on external data (the NSIDC ice front netCDF used to
% derive catchment boundaries, exactly as in runme_greenibis.m).

%% CONFIGURATION - edit for your environment {{{
domain          = 'CW';
fric            = 1;

% Catchment ids to analyze
if strcmp(domain,'CW')
	catchment_ids   = [3 4 5 6 7 8 9 10 72 81 82 85 218 219 220 221 222];
elseif strcmp(domain,'NW')
   catchment_ids = [1 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 43 44
5 46 47 48 49 78 86 103 104 105 127 147 148 160 162 164 165 167 172 174 175 178
180 181 183 185 195 196 197 207 208 214 225 238 239 240 242];
end

% Model paths (before softening / after softening, not re-solved / after
% softening, re-solved) - same repository/prefix/step naming as
% runme_snapshot.m's organizer, so these follow `domain` automatically.
model_folder         = './Models_snapshot';
model_prefix         = 'Model_';
model_before_path    = fullfile(model_folder, [model_prefix domain '_Friction.mat']);
model_after_path  = fullfile(model_folder, [model_prefix domain '_Rheology.mat']);

% Same external data used by runme_snapshot.m to build catchment boundaries
nsidc_icefronts_netcdf = '/home/badgeley/ModelData/Greenland/IceFrontsGreene/NSIDC-0793_19720915-20220215_V01.0.nc';
exp_folder              = './Exp';
contour_resolution_m    = 200; % simplification tolerance, same as runme_snapshot.m

% Same thresholds used in runme_snapshot.m's shear margin softening loop
min_thickness_m  = 10;
min_velocity_myr = 50;
rsquare_threshold = 0.8;
softening_vel_percentile = 95;

% Where to write the per-catchment .mat exports for the Python script
output_dir = './shear_margin_softening_exports';
%}}}

%% SETUP {{{
addpath(fullfile(getenv('ISSM_DIR'), 'src', 'm', 'dev'));
devpath;

if ~exist(exp_folder, 'dir')
    mkdir(exp_folder);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% output_dir is a relative path, so it resolves against MATLAB's current
% directory - print the resolved absolute path so it's obvious where
% exports are actually landing (e.g. if this is run from the wrong cwd,
% it's easy to end up silently nesting a new output_dir inside an old one).
fprintf('Exporting to: %s\n', fullfile(pwd, output_dir));

disp('Loading before/rheology/after models');
model_before   = loadmodel(model_before_path);
model_after    = loadmodel(model_after_path);

nv = model_before.mesh.numberofvertices;
if model_after.mesh.numberofvertices ~= nv
    error(['model_before and model_after must share the same mesh ' ...
        '(numberofvertices differs) - are these really the before/after ' ...
        'results for the same domain/run?']);
end

vel_before_full = model_before.results.StressbalanceSolution.Vel;
vel_after_full  = model_after.results.StressbalanceSolution.Vel;
rheology_before_full = model_before.materials.rheology_B;
rheology_after_full  = model_after.materials.rheology_B;
%}}}

%% PER-CATCHMENT ANALYSIS {{{
for catchment = catchment_ids

    fprintf('\n== Catchment %d ==\n', catchment);

    %-- Get (or build) the catchment boundary, same method as runme_greenibis.m {{{
    exp_file = fullfile(exp_folder, ['MouginotGreene_extended_catchment_' num2str(catchment) '.exp']);
    if ~isfile(exp_file)
        disp('   -- Building catchment boundary from NSIDC ice front catchments');
        X = ncread(nsidc_icefronts_netcdf, 'x');
        Y = ncread(nsidc_icefronts_netcdf, 'y');
        catchmentID = double(ncread(nsidc_icefronts_netcdf, 'catchment_name'));

        MASK = double(catchmentID == catchment);
        c = contourc(X, Y, MASK', [0.5 0.5]);

        done = 0; ii = 1; jj = 1; s = [];
        while (ii < length(c))
            num = c(2, ii); ii = ii + 1;
            temp = [c(1, ii:(ii + num - 1))' c(2, ii:(ii + num - 1))'];
            temp = dpsimplify(temp, contour_resolution_m);
            s(jj).x = temp(:, 1); s(jj).y = temp(:, 2);
            ii = ii + num; jj = jj + 1;
        end
        expwrite(s, exp_file);
    end

    flag_nodes = ContourToMesh(model_before.mesh.elements, model_before.mesh.x, model_before.mesh.y, exp_file, 'node', 1);
    pos_flag_nodes = find(flag_nodes == 1);
    fprintf('   %d flagged nodes found for basin %d\n', length(pos_flag_nodes), catchment);
    if isempty(pos_flag_nodes)
        warning('No nodes found for catchment %d - skipping.', catchment);
        continue
    end
    flag_elements = sum(ismember(model_before.mesh.elements, pos_flag_nodes), 2);
    flag_elements(flag_elements > 0) = 1;
    %}}}

    %-- Extract the catchment sub-model (same as runme_greenibis.m) {{{
    mde = extract(model_before, flag_elements, 'spccheck', 0);
    mde = mechanicalproperties(mde, mde.inversion.vx_obs, mde.inversion.vy_obs);
    %}}}

    %-- Strain rate, velocity error, and velocity error slope (same as runme_snapshot.m) {{{
    ice_levelset = mde.mask.ice_levelset;

    analysis_mask = ice_levelset;
    analysis_mask(analysis_mask > 0) = NaN;
    analysis_mask(analysis_mask < 0) = 1;
    analysis_mask(mde.geometry.thickness < min_thickness_m) = NaN;
    analysis_mask(mde.results.StressbalanceSolution.Vel < min_velocity_myr) = NaN;

    v_error = mde.results.StressbalanceSolution.Vel - mde.initialization.vel;
    [~, ~, v_error_slope_elements] = slope(mde, v_error);
    v_error_slope = averaging(mde, abs(v_error_slope_elements), 0);
    strain_rate = averaging(mde, mde.results.strainrate.effectivevalue, 0);

    % Velocity cutoff used below to exclude the fastest ice from softening -
    % an unweighted percentile over vertices, matching the actual
    % runme_greenibis.m pipeline that produced these models (not an
    % area-weighted percentile - the mesh is much denser at the shear
    % margins than in the catchment interior, so this cutoff skews toward
    % the margins; that's called out explicitly in the summary PDF).
    vel_cutoff = prctile(mde.results.StressbalanceSolution.Vel, softening_vel_percentile);
    fprintf('   shear margin maximum velocity cutoff (unweighted %gth percentile): %g\n', softening_vel_percentile, vel_cutoff);
    %}}}

    %-- Linear fit of strain rate vs. velocity error slope (same relationship as runme_snapshot.m) {{{
    % NOTE: runme_snapshot.m calls `[fitobj,gof,~]=polyfit(X,Y,1)` and reads
    % gof.rsquared, which is not a field base MATLAB's polyfit returns - that
    % call is only valid with the Curve Fitting Toolbox's `fit`/gof pattern
    % (as used in strainrate_vs_error.m). Since fenris has no toolboxes
    % installed, this uses base MATLAB's polyfit/polyval instead, with R^2
    % and RMSE computed by hand (see local function fitstats_poly1 below).
    pos_fit = find(~isnan(v_error_slope .* analysis_mask) & ~isnan(strain_rate .* analysis_mask));
    X = v_error_slope(pos_fit) .* analysis_mask(pos_fit);
    Y = strain_rate(pos_fit) .* analysis_mask(pos_fit);

    valid_fit = false(mde.mesh.numberofvertices, 1);
    valid_fit(pos_fit) = true;
    softened = false(mde.mesh.numberofvertices, 1);

    if length(pos_fit) < 3
        warning('Catchment %d has too few valid points (%d) for a fit - skipping softening classification.', catchment, length(pos_fit));
        fit_slope = NaN; fit_intercept = NaN; fit_rsquare = NaN; fit_rmse = NaN; fit_n = length(pos_fit);
        was_softened = false;
    else
        [fit_slope, fit_intercept, fit_rsquare, fit_rmse] = fitstats_poly1(X, Y);
        fit_n = length(X);

        fprintf('   r-squared of strain rate vs. velocity error slope fit: %g\n', fit_rsquare);
        fprintf('   slope: %g, intercept: %g\n', fit_slope, fit_intercept);

        was_softened = fit_rsquare >= rsquare_threshold;
        if ~was_softened
            disp('   r-squared below threshold - shear margins were not softened for this catchment.');
        else
            strain_rate_fit = (v_error_slope .* fit_slope) + fit_intercept;
            if fit_slope >= 0
                pos1 = find(strain_rate >= strain_rate_fit);
            else
                pos1 = find(strain_rate <= strain_rate_fit);
            end
            pos2 = find(mde.results.StressbalanceSolution.Vel < vel_cutoff);

            pos = intersect(intersect(pos1, pos2), pos_fit);
            softened(pos) = true;
        end
    end
    %}}}

    %-- Before/after rheology and velocity, mapped back through extractedvertices {{{
    gv = mde.mesh.extractedvertices; % local -> global vertex index map

    rheology_before = rheology_before_full(gv);
    rheology_after  = rheology_after_full(gv);
    rheology_diff   = rheology_after - rheology_before;

    vel_before = vel_before_full(gv);
    vel_after  = vel_after_full(gv);
    vel_diff   = vel_after - vel_before;
    %}}}

    %-- Export plain arrays for Python {{{
    x = mde.mesh.x;
    y = mde.mesh.y;
    elements = mde.mesh.elements; % 1-based triangle connectivity

    out_path = fullfile(output_dir, ['catchment_' num2str(catchment) '.mat']);
    save(out_path, 'catchment', 'x', 'y', 'elements', 'ice_levelset', ...
        'strain_rate', 'v_error', 'v_error_slope', ...
        'rheology_before', 'rheology_after', 'rheology_diff', ...
        'vel_before', 'vel_after', 'vel_diff', ...
        'valid_fit', 'softened', 'was_softened', ...
        'fit_slope', 'fit_intercept', 'fit_rsquare', 'fit_rmse', 'fit_n', ...
        'min_thickness_m', 'min_velocity_myr', 'rsquare_threshold', ...
        'vel_cutoff', 'softening_vel_percentile', '-v7');
    fprintf('   Wrote %s\n', out_path);
    %}}}
end
%}}}

disp('Done.');

%% Local functions {{{
function [slope_, intercept_, rsquare_, rmse_] = fitstats_poly1(x, y) % {{{
    %FITSTATS_POLY1 - degree-1 polynomial fit with R^2 and RMSE, using only
    %base MATLAB (polyfit/polyval), as a toolbox-free replacement for the
    %Curve Fitting Toolbox's fit(x,y,'poly1')/gof pattern.
    x = x(:); y = y(:);
    p = polyfit(x, y, 1);
    slope_ = p(1);
    intercept_ = p(2);

    yhat = polyval(p, x);
    residuals = y - yhat;
    ss_res = sum(residuals.^2);
    ss_tot = sum((y - mean(y)).^2);
    rsquare_ = 1 - ss_res / ss_tot;

    n = length(x);
    dof = max(n - numel(p), 1); % degrees of freedom (n - number of fit coefficients)
    rmse_ = sqrt(ss_res / dof);
end % }}}
%}}}
