function [shear, shear_fraction, effective] = lateral_shear_strain_rate(md, vx, vy)
%LATERAL_SHEAR_STRAIN_RATE: Rotates the horizontal strain rate tensor of a
%velocity field into flow-aligned coordinates and returns the lateral
%(across-flow) shear component on the vertices.
%
%   The effective strain rate mixes lateral shear with along-flow stretching
%   and across-flow convergence, so it is high in the trunk near a terminus as
%   well as in the shear margins. The lateral shear component isolates the
%   margins: with s along flow and n across flow,
%       e_sn = (eyy - exx)*sin*cos + exy*(cos^2 - sin^2)
%   which is half the across-flow gradient of the along-flow speed (plus a
%   flow-curvature term). It changes sign across the centerline, so the signed
%   value is averaged onto the vertices before taking the magnitude - that
%   keeps the trunk center at ~0 instead of smearing both margins into it.
%
%   INPUTS:
%           md: ISSM model (2d mesh)
%           vx, vy: velocity components (m/yr) on the vertices, e.g.
%               md.inversion.vx_obs and md.inversion.vy_obs
%   OUTPUTS:
%           shear: |lateral shear strain rate| (1/yr) on the vertices
%           shear_fraction: shear ./ effective, from 0 (no lateral shear) to
%               1 (pure lateral shear)
%           effective: effective strain rate (1/yr) on the vertices, same
%               definition as mechanicalproperties' strainrate.effectivevalue
%   EXAMPLE:
%           [shear, shear_fraction] = lateral_shear_strain_rate(md, ...
%              md.inversion.vx_obs, md.inversion.vy_obs);
md = mechanicalproperties(md, vx, vy);
exx = md.results.strainrate.xx;
eyy = md.results.strainrate.yy;
exy = md.results.strainrate.xy;

% Flow direction on each element (undefined where the ice is not moving)
vx_elements = mean(vx(md.mesh.elements), 2);
vy_elements = mean(vy(md.mesh.elements), 2);
speed_elements = hypot(vx_elements, vy_elements);
c = vx_elements ./ speed_elements;
s = vy_elements ./ speed_elements;
c(speed_elements == 0) = 1;
s(speed_elements == 0) = 0;

shear_elements = (eyy - exx) .* s .* c + exy .* (c.^2 - s.^2);

shear = abs(averaging(md, shear_elements, 0));
effective = averaging(md, md.results.strainrate.effectivevalue, 0);

shear_fraction = shear ./ effective;
shear_fraction(effective == 0) = 0;
end
