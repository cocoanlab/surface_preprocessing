function [flist, search_path, error_flag] = search_files(fname, maxdepth)

% The function search files from the specified path and name.
%
%
% :Usage:
% ::
%     [flist, search_path, error_flag] = search_files(fname, maxdepth)
%
%
% :Input:
% ::
%   - fname              search path to files.
%   - maxdepth           maximum depth of search.
%
%
% :Output:
% ::
%     flist              list of searched files.
%     search_path        actual search path under maximal depth.
%     error_flag         success/fail of search (0: success, 1: fail).
%
%
% :Example:
% ::
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Jan 2019  Choong-Wan Woo & Jae-Joong Lee
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..


error_flag = 0;
flist = [];

[pathstr, name, ext] = fileparts(fname);
for depth_i = 1:maxdepth
    search_path = [pathstr '/' repmat('*/', 1, depth_i - 1) name ext];
    dirout = dir(search_path);
    if ~isempty(dirout)
        break
    else
        if depth_i == maxdepth
            warning('Failed to search files.');
            error_flag = 1;
            return;
        end
    end
end
flist = cellstr(strcat(cat(1, dirout.folder), '/', cat(1, dirout.name)));

end