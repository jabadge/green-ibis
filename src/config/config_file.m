function [p, data] = config_file()
%CONFIG_FILE  Location and contents of the shared path config file.
%
%   The same file the Python side reads and writes (src/config/config.py), so a
%   path answered once is known to both languages. Location:
%      $GREENIBIS_CONFIG                        if set, else
%      $XDG_CONFIG_HOME/green-ibis/paths.json    if set, else
%      ~/.config/green-ibis/paths.json
%
%   Usage:
%      p = config_file();
%      [p, data] = config_file();
%
%   Output:
%      p    : (char)   full path to the config file (may not exist yet)
%      data : (struct) its decoded contents (empty struct if absent/unreadable)

override = getenv('GREENIBIS_CONFIG');
if ~isempty(override)
    p = char(override);
else
    xdg = getenv('XDG_CONFIG_HOME');
    if isempty(xdg)
        xdg = fullfile(home_dir(), '.config');
    end
    p = fullfile(char(xdg), 'green-ibis', 'paths.json');
end

data = struct();
if nargout > 1 && isfile(p)
    try
        data = jsondecode(fileread(p));
    catch
        data = struct();   % unreadable or malformed: treat as empty
    end
end

end
