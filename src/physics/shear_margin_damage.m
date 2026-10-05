function [damage, margin, info, shear_peak] = shear_margin_damage(x, y, elements, vx, vy, shear, shear_fraction, valid, varargin)
%SHEAR_MARGIN_DAMAGE: Finds the shear margins of a glacier from its lateral
%shear strain rate and returns a softening factor ("damage") for them.
%
%   A margin is defined relative to the LOCAL peak of the lateral shear, taken
%   across flow: for every vertex, shear_peak is the highest shear found within
%   'peak_radius' to either side of it (across flow) and 'peak_along' up or
%   down flow. A vertex belongs to a margin when its shear is a large enough
%   share of that peak, so the edge of a margin follows the shape of the shear
%   profile at that position along the glacier - its width does not grow and
%   shrink with how strong the shear happens to be there.
%
%   The margins are found with two thresholds joined through the mesh
%   connectivity:
%     - seeds:  shear-dominated vertices close to the local peak
%     - grown:  shear-dominated vertices further down the flanks of the peak
%   A margin is a connected patch (along mesh edges) of grown vertices that
%   contains at least one seed and covers at least 'min_area'.
%
%   Every vertex has a local peak, even in ice that is barely shearing, so
%   margins are only looked for where the shear is at least
%   'damage_shear_min' - the shear at which the softening starts. Nothing here
%   depends on the extent of the mesh, so it can be run over a whole domain at
%   once: every glacier is treated the same way, and ice shared by
%   neighboring catchments is only softened once.
%
%   The softening is set from the absolute lateral shear, not the relative
%   one, as an enhancement factor E on the rate factor A of Glen's flow law:
%   E grows with the shear (as a power law) from 1 at 'damage_shear_min' to
%   'max_enhancement' at 'damage_shear_max', and is 1 outside the margins.
%   Slow glaciers therefore get little or no softening without a separate
%   test.
%
%   ISSM's rheology_B is A^(-1/n), so multiplying A by E divides rheology_B
%   by E^(1/n) - that divisor is the damage returned here (10^(1/3) ~ 2.15
%   at most with the defaults).
%
%   NOTE on the amount of softening: the default maximum of E = 10 comes from
%   Sun et al. (2017), "Ice shelf fracture parameterization in an ice sheet
%   model", The Cryosphere, 11, 2543-2554, doi:10.5194/tc-11-2543-2017, who
%   found that heavily damaged ice (damage ~ 1/2, equivalent to E = 8) can be
%   emulated by an order of magnitude change in A (E ~ 10). Caveats:
%     - That result is for crevasse damage on an IDEALIZED ICE SHELF
%       (MISMIP+), where the modeled damage on grounded ice was small.
%       Applying it to the shear margins of grounded ice is an extrapolation.
%     - The paper only motivates the maximum. The power law growth of E with
%       the lateral shear strain rate, and the strain rates at which it
%       starts and saturates ('damage_shear_min'/'damage_shear_max'), are
%       choices made here. Slower glaciers never reach the full factor.
%
%   INPUTS:
%           x, y: vertex coordinates (m)
%           elements: triangle connectivity (1-based)
%           vx, vy: velocity components on the vertices, used only for the
%               flow direction (same field the shear was computed from)
%           shear: |lateral shear strain rate| (1/yr) on the vertices, see
%               lateral_shear_strain_rate
%           shear_fraction: lateral shear ./ effective strain rate (0 to 1)
%           valid: logical, vertices that may be part of a margin (ice
%               covered, thick enough, fast enough)
%   NAME/VALUE OPTIONS:
%           'seed_shear_fraction': minimum shear fraction of a seed (0.8)
%           'grow_shear_fraction': minimum shear fraction of a grown vertex (0.6)
%           'seed_peak_ratio': minimum shear ./ shear_peak of a seed (0.8)
%           'grow_peak_ratio': minimum shear ./ shear_peak of a grown vertex,
%               i.e. where the edge of a margin sits on the flank of the
%               peak (0.5)
%           'peak_radius': how far across flow to look for the peak (m) (3000)
%           'peak_along': how far along flow to look for the peak (m) (1000)
%           'min_area': smallest margin patch to keep (m^2) (2e6)
%           'max_enhancement': largest enhancement factor on A (10)
%           'rheology_n': Glen's flow law exponent, a scalar or one value
%               per vertex (3)
%           'damage_shear_min': lateral shear (1/yr) where softening starts,
%               and the least shear anywhere in a margin (0.03)
%           'damage_shear_max': lateral shear (1/yr) where the enhancement
%               reaches max_enhancement (1)
%   OUTPUTS:
%           damage: divisor of rheology_B on the vertices (>= 1), equal to
%               enhancement.^(1/n); the softened rheology is
%               rheology_B ./ damage
%           margin: 0 outside the margins, 1 grown, 2 seed
%           info: struct with n_margins, max_damage (the largest possible
%               divisor of rheology_B), and the options used
%           shear_peak: local across-flow peak of the lateral shear (1/yr) on
%               the vertices (0 where not valid)
%   EXAMPLE:
%           [shear, shear_fraction] = lateral_shear_strain_rate(md, ...
%              md.inversion.vx_obs, md.inversion.vy_obs);
%           valid = md.mask.ice_levelset < 0 & md.geometry.thickness >= 10 ...
%              & md.inversion.vel_obs >= 50;
%           damage = shear_margin_damage(md.mesh.x, md.mesh.y, ...
%              md.mesh.elements, md.inversion.vx_obs, md.inversion.vy_obs, ...
%              shear, shear_fraction, valid);
%           md.materials.rheology_B = md.materials.rheology_B ./ damage;
p = inputParser;
p.addParameter('seed_shear_fraction', 0.8);
p.addParameter('grow_shear_fraction', 0.6);
p.addParameter('seed_peak_ratio', 0.8);
p.addParameter('grow_peak_ratio', 0.5);
p.addParameter('peak_radius', 3000);
p.addParameter('peak_along', 1000);
p.addParameter('min_area', 2e6);
p.addParameter('max_enhancement', 10);
p.addParameter('rheology_n', 3);
p.addParameter('damage_shear_min', 0.03);
p.addParameter('damage_shear_max', 1);
p.parse(varargin{:});
opt = p.Results;

nv = numel(x);
x = x(:); y = y(:); vx = vx(:); vy = vy(:); shear = shear(:); shear_fraction = shear_fraction(:);
valid = logical(valid(:)) & isfinite(shear) & isfinite(shear_fraction);

% Area belonging to each vertex (a third of each surrounding triangle)
x_el = x(elements); y_el = y(elements);
element_area = 0.5 * abs((x_el(:,2) - x_el(:,1)) .* (y_el(:,3) - y_el(:,1)) ...
    - (y_el(:,2) - y_el(:,1)) .* (x_el(:,3) - x_el(:,1)));
vertex_area = accumarray(elements(:), repmat(element_area / 3, 3, 1), [nv 1]);

damage = ones(nv, 1);
margin = zeros(nv, 1);
info = opt;
info.rheology_n = mean(opt.rheology_n(:));
info.max_damage = opt.max_enhancement ^ (1 / min(opt.rheology_n(:)));
info.n_margins = 0;
shear_peak = zeros(nv, 1);
if ~any(valid)
    return
end

% Local peak of the shear across flow. For each valid vertex, look at the
% valid vertices inside a box centered on it that reaches peak_radius across
% flow and peak_along along flow (done in chunks to keep the distance
% matrices small).
pos = find(valid);
xv = x(pos); yv = y(pos); shear_v = shear(pos);
speed = hypot(vx(pos), vy(pos));
c = vx(pos) ./ speed; s = vy(pos) ./ speed;
c(speed == 0) = 1; s(speed == 0) = 0;
peak_v = zeros(numel(pos), 1);
chunk = 500;
for first = 1:chunk:numel(pos)
    k = (first:min(first + chunk - 1, numel(pos)))';
    dx = xv' - xv(k);
    dy = yv' - yv(k);
    in_box = abs(dx .* c(k) + dy .* s(k)) <= opt.peak_along ...
        & abs(-dx .* s(k) + dy .* c(k)) <= opt.peak_radius;
    peak_v(k) = max(in_box .* shear_v', [], 2);
end
shear_peak(pos) = peak_v;

strong = valid & shear >= opt.damage_shear_min;
seed = strong & shear_fraction >= opt.seed_shear_fraction & shear >= opt.seed_peak_ratio * shear_peak;
grown = strong & shear_fraction >= opt.grow_shear_fraction & shear >= opt.grow_peak_ratio * shear_peak;

% Connected patches of grown vertices, joined along mesh edges
edges = [elements(:, [1 2]); elements(:, [2 3]); elements(:, [3 1])];
edges = edges(grown(edges(:,1)) & grown(edges(:,2)), :);
patch_id = conncomp(graph(edges(:,1), edges(:,2), [], nv))';

patch_area = accumarray(patch_id(grown), vertex_area(grown), [max(patch_id) 1]);
keep = unique(patch_id(seed));
keep = keep(patch_area(keep) >= opt.min_area);
in_margin = grown & ismember(patch_id, keep);

margin(in_margin) = 1;
margin(in_margin & seed) = 2;

n = opt.rheology_n(:) .* ones(nv, 1);
ramp = log(shear(in_margin) / opt.damage_shear_min) / log(opt.damage_shear_max / opt.damage_shear_min);
damage(in_margin) = opt.max_enhancement .^ (min(max(ramp, 0), 1) ./ n(in_margin));

info.n_margins = numel(keep);
end
