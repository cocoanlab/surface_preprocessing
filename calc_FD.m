function FD = calc_FD(mvmt_param, param_software)

% The function calculates framewise displacement (FD) based on Power et al., 2012, Neuroimage.
%
%
% :Usage:
% ::
%     FD = calc_FD(mvmt_param, varargin)
%
%
% :Input:
% ::
%   - mvmt_param         the movement parameter (t x 6) without zero-padding
%                        (it will be done in this function).
%   - param_software     the software that was used for calculating
%                        the movement parameter.
%                        (e.g., 'AFNI', 'FSL', 'SPM')
%
%
% :Optional Input:
% ::
%
%
% :Output:
% ::
%     FD
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


radius_brain = 50; % Based on Power's assumption

if size(mvmt_param, 2) ~= 6
    error('The number of movement parameter variables should be 6.');
end

switch param_software
    case 'AFNI' % the order of translation and rotation is reversed!
        mvmt_trans = mvmt_param(:,4:6);
        mvmt_rot = mvmt_param(:,1:3) * pi / 180; % degree
    case 'FSL'
        mvmt_trans = mvmt_param(:,4:6);
        mvmt_rot = mvmt_param(:,1:3); % radian
    case 'SPM'
        mvmt_trans = mvmt_param(:,1:3);
        mvmt_rot = mvmt_param(:,4:6); % radian
end
mvmt_rot = mvmt_rot * radius_brain;
mvmt_all = [mvmt_trans, mvmt_rot];
mvmt_all_diff = [zeros(1, size(mvmt_all, 2)); diff(mvmt_all)];
FD = sum(abs(mvmt_all_diff), 2);

end