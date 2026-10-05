%SHEAR_MARGIN_SOFTENING_ANALYSIS
%
% Finds the shear margins of a whole domain from the lateral shear strain
% rate of the OBSERVED velocities (see lateral_shear_strain_rate and
% shear_margin_damage), and exports them, with their softening factor
% ("damage"), as a plain (non-classdef) .mat file for
% plot_shear_margin_summary.py to map.
%
% This is the same definition the 'Rheology' step of runme_snapshot.m applies:
%
%   1. rotate the strain rate tensor into flow-aligned coordinates and keep
%      the lateral shear component (the trunk and the stretching ice near the
%      terminus have little of it, the margins are dominated by it)
%   2. find the local peak of that shear across flow, seed the margins close
%      to the peak, then grow them through the mesh connectivity down its
%      flanks, so the margin edge follows the local shape of the shear profile
%   3. scale an enhancement factor on the rate factor A with the lateral
%      shear strain rate itself (up to ~10, after Sun et al., 2017), and
%      convert it to the matching divisor of rheology_B (the "damage",
%      at most 10^(1/3) ~ 2.15). See the NOTE in shear_margin_damage for
%      the reference and for the caveats: the factor of 10 is taken from an
%      idealized ice shelf study and applied here to grounded shear margins,
%      and the power law growth with shear is a choice made here.
%
% This script does not modify or re-save a model. It works from:
%
%   model_before    - the model BEFORE shear margin softening
%                      (runme_snapshot.m 'Friction' step output)
%   model_after     - optional: the same mesh AFTER the softening of the
%                      runme_snapshot.m 'Rheology' step. If the file exists,
%                      the softening actually in it (rheology_B before /
%                      after) is exported too, to compare against the damage
%                      computed here. They differ only by the smoothing that
%                      step applies to rheology_B afterwards.
%
% Usage:
%   Run setup_paths, edit the CONFIGURATION block below to match your
%   environment, then:
%       shear_margin_softening_analysis

%% CONFIGURATION - edit for your environment {{{
domain          = 'CW';

% Model paths - same repository/prefix/step naming as runme_snapshot.m's
% organizer, so these follow `domain` automatically.
model_folder         = './Models_snapshot';
model_prefix         = 'Model_';
model_before_path    = fullfile(model_folder, [model_prefix domain '_Friction.mat']);
model_after_path     = fullfile(model_folder, [model_prefix domain '_Rheology.mat']);

% Vertices that may be part of a shear margin: ice covered, at least this
% thick, and with an OBSERVED speed of at least this much
min_thickness_m  = 10;
min_velocity_myr = 50;

% Shear margin definition and damage scaling - see shear_margin_damage for
% what each one does
margin_options = { ...
    'seed_shear_fraction', 0.8, ...   % lateral shear / effective strain rate of a seed
    'grow_shear_fraction', 0.6, ...   % ... and of the ice the margin grows into
    'seed_peak_ratio',     0.8, ...   % lateral shear / local across-flow peak shear of a seed
    'grow_peak_ratio',     0.5, ...   % ... and of the ice the margin grows into (the margin edge)
    'peak_radius',         3000, ...  % how far across flow to look for the local peak (m)
    'peak_along',          1000, ...  % how far along flow to look for the local peak (m)
    'min_area',            2e6, ...   % smallest margin patch to keep (m^2)
    'max_enhancement',     10, ...    % largest enhancement factor on the rate factor A (Sun et al., 2017)
    'damage_shear_min',    0.03, ...  % lateral shear (1/yr) where softening starts; no margins below it
    'damage_shear_max',    1};        % lateral shear (1/yr) where it reaches max_enhancement

% Where to write the export for the Python script
output_file = ['./shear_margin_softening_' domain '.mat'];
%}}}

%% SETUP {{{
addpath(fullfile(getenv('ISSM_DIR'), 'src', 'm', 'dev'));
devpath;

disp('Loading the model before softening');
model_before = loadmodel(model_before_path);
%}}}

%% SHEAR MARGINS AND DAMAGE {{{
disp('Finding the shear margins of the whole domain');

%-- Lateral shear of the observed velocities
vel_obs = hypot(model_before.inversion.vx_obs, model_before.inversion.vy_obs);
[lateral_shear, shear_fraction] = lateral_shear_strain_rate(model_before, ...
    model_before.inversion.vx_obs, model_before.inversion.vy_obs);

valid = model_before.mask.ice_levelset < 0 & model_before.geometry.thickness >= min_thickness_m ...
    & vel_obs >= min_velocity_myr;

%-- Shear margins and their damage
[damage, margin, margin_info] = shear_margin_damage(model_before.mesh.x, model_before.mesh.y, ...
    model_before.mesh.elements, model_before.inversion.vx_obs, model_before.inversion.vy_obs, ...
    lateral_shear, shear_fraction, valid, margin_options{:}, ...
    'rheology_n', averaging(model_before, model_before.materials.rheology_n, 0));

fprintf('   %d margin patches, %d vertices, %d softened (damage > 1), mean damage %g\n', ...
    margin_info.n_margins, nnz(margin), nnz(damage > 1), mean(damage(margin > 0)));
%}}}

%% SOFTENING ALREADY IN THE SAVED MODEL, IF THERE IS ONE {{{
damage_applied = [];
if isfile(model_after_path)
    disp('Loading the model after softening');
    model_after = loadmodel(model_after_path);
    if model_after.mesh.numberofvertices ~= model_before.mesh.numberofvertices
        error(['model_before and model_after must share the same mesh ' ...
            '(numberofvertices differs) - are these really the before/after ' ...
            'results for the same domain/run?']);
    end
    damage_applied = model_before.materials.rheology_B ./ model_after.materials.rheology_B;
end
%}}}

%% EXPORT {{{
x = model_before.mesh.x;
y = model_before.mesh.y;
elements = model_before.mesh.elements; % 1-based triangle connectivity
ice_levelset = model_before.mask.ice_levelset;

save(output_file, 'domain', 'x', 'y', 'elements', 'ice_levelset', ...
    'margin', 'damage', 'damage_applied', 'margin_info', '-v7');
fprintf('Wrote %s\n', fullfile(pwd, output_file));
%}}}
