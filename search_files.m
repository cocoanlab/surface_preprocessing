function [flist, search_path, error_flag] = search_files(fname, varargin)

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
%
%
% :Optional Input:
% ::
%   - maxdepth           maximum depth of search. (default: 1)
%   - char               output char array. (default: cell array)
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
maxdepth = 1;
dochar = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'maxdepth'}
                maxdepth = varargin{i+1};
            case {'char'}
                dochar = true;
        end
    end
end

    
flist = [];
[pathstr, name, ext] = fileparts(fname);
for depth_i = 1:maxdepth
    search_path = [pathstr '/' repmat('*/', 1, depth_i - 1) name ext];
    [~, flist_part] = system(['for file in ' search_path '; do echo $file; done']);
    flist_part(end) = []; % delete the last '\n'
    if strcmp(flist, search_path) % cannot find files, so only the search path was out
        if depth_i == maxdepth
            warning('Failed to search files.');
            error_flag = 1;
            flist = [];
            return;
        end
    end
    flist = [flist flist_part];
end

flist = strsplit(flist, '\n');
flist = flist(:); % vectorize
if dochar
    flist = char(flist);
end

end