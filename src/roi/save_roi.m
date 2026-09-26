function roi = save_roi(roi, varargin)
%SAVE_ROI  Save an ROI struct to rois/<name>.mat and an ISSM .exp contour.
%
%   Usage:
%      roi = save_roi(roi)
%      roi = save_roi(roi, 'dir', '/some/folder')
%
%   Input:
%      roi : ROI struct with .name, .x, .y (see select_roi)
%   Name/value options:
%      'dir' : (char) output folder (default: project_config().roi_dir)
%
%   Output:
%      roi : same struct with .file set to the saved .mat path

p = inputParser;
p.addParameter('dir','', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
odir = char(p.Results.dir);
if isempty(odir)
    odir = data_path('roi_dir');
end
if ~isfolder(odir), mkdir(odir); end

if ~isfield(roi,'name') || isempty(roi.name)
    roi.name = ['roi_' datestr(now,'yyyymmdd_HHMMSS')];
end

matfile = fullfile(odir, [roi.name '.mat']);
roi.file = matfile;
save(matfile, 'roi');

% Also write an ISSM Argus .exp contour (usable by ContourToNodes / SectionValues / meshing)
expfile = fullfile(odir, [roi.name '.exp']);
try
    c = struct('x', roi.x(:), 'y', roi.y(:), 'density', 1);
    expwrite(c, expfile);
catch ME
    warning('save_roi:exp', 'Could not write .exp (%s): %s', expfile, ME.message);
end

fprintf('Saved ROI "%s" -> %s (+ .exp)\n', roi.name, matfile);

end
