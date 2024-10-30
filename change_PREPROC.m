function change_PREPROC(PREPROC_path, orig_basedir, new_basedir, searchdepth)

% The function change basedir of PREPROC.
%
%
% :Usage:
% ::
%     change_PREPROC(PREPROC_path, orig_basedir, new_basedir)
%
%
% :Input:
% ::
%   - PREPROC_path       search path to PREPROC.mat
%   - orig_basedir       original basedir of PREPROC.mat
%   - new_basedir        new basedir of PREPROC.mat
%
%
% :Optional Input:
% ::
%
%
% :Output:
% ::
%
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


if nargin < 4; searchdepth = 1; end
fprintf('\n');
PREPROC_list = search_files(PREPROC_path, 'maxdepth', searchdepth);
for i = 1:numel(PREPROC_list)
    fprintf('Converting basedir of PREPROC.mat:   %s\n', PREPROC_list{i});
    load(PREPROC_list{i}, 'PREPROC');
    PREPROC = replaceSubstrInUnknownDepth(PREPROC, orig_basedir, new_basedir);
    save(PREPROC_list{i}, 'PREPROC');
end
fprintf('\n');

end

function dataOut = replaceSubstrInUnknownDepth(dataIn, oldStr, newStr)
    % Check the type of the current element
    if ischar(dataIn)
        % If it's a string, perform the replacement
        dataOut = strrep(dataIn, oldStr, newStr);
    elseif isstruct(dataIn)
        % If it's a struct array, recursively apply to each element
        if numel(dataIn) > 1
            dataOut = dataIn;  % Initialize output struct array
            for i = 1:numel(dataIn)
                dataOut(i) = replaceSubstrInUnknownDepth(dataIn(i), oldStr, newStr);
            end
        else
            % If it's a single struct, process each field
            fields = fieldnames(dataIn);
            dataOut = struct();  % Initialize output struct
            for i = 1:numel(fields)
                field = fields{i};
                dataOut.(field) = replaceSubstrInUnknownDepth(dataIn.(field), oldStr, newStr);
            end
        end
    elseif iscell(dataIn)
        % If it's a cell array, recursively apply the function to each element
        dataOut = cell(size(dataIn)); % Initialize the output cell array
        for i = 1:numel(dataIn)
            dataOut{i} = replaceSubstrInUnknownDepth(dataIn{i}, oldStr, newStr);
        end
    else
        % For other types (e.g., numbers), leave unchanged
        dataOut = dataIn;
    end
end