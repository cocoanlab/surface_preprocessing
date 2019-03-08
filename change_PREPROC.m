function change_PREPROC(PREPROC_path, orig_basedir, new_basedir)

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


fprintf('\n');
PREPROC_list = search_files(PREPROC_path, 'maxdepth', 3);
for i = 1:numel(PREPROC_list)
    fprintf('Converting basedir of PREPROC.mat:   %s\n', PREPROC_list{i});
    load(PREPROC_list{i}, 'PREPROC');
    dat = struct2cell(PREPROC);
    wh_numeric = cellfun(@isnumeric, dat);
    wh_struct = cellfun(@isstruct, dat);
    dat(~wh_struct & ~wh_numeric) = cellfun(@(x) strrep(x, orig_basedir, new_basedir), dat(~wh_struct & ~wh_numeric), 'UniformOutput', false);
    dat(wh_struct) = cellfun(@(y) structfun(@(x) strrep(x, orig_basedir, new_basedir), y, 'UniformOutput', false), dat(wh_struct), 'UniformOutput', false);
    PREPROC = cell2struct(dat, fieldnames(PREPROC));
    save(PREPROC_list{i}, 'PREPROC');
end
fprintf('\n');

end