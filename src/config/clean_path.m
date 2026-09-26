function p = clean_path(s)
%CLEAN_PATH  Tidy a pasted or drag-and-dropped path.
%
%   Strips surrounding whitespace and quotes (Finder/Explorer drag-and-drop
%   often wraps a path in them) and expands a leading ~ to the home folder.
%
%   Usage:
%      p = clean_path('  "~/data/velocity"  ')
%
%   Input:
%      s : (char/string) the raw text
%   Output:
%      p : (char) the tidied path

p = char(s);
p = strtrim(p);
p = regexprep(p, '^["'']|["'']$', '');
p = strtrim(p);
if startsWith(p, '~')
    p = fullfile(home_dir(), p(2:end));
end
end
