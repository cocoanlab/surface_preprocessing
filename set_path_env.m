function [rootdir, basedir, rscdir, gitdir] = set_path_env(varargin)

% The function gets useful directory names and sets basic environment
% (path to toolboxes, environment variables, maximum number of CPU cores)
% depending on current hostname.
%
%
% :Usage:
% ::
%     [rootdir, basedir, rscdir, gitdir] = set_path_env(varargin)
%
%
% :Optional Input:
% ::
%   - hostname           custom hostname.
%                        if this option is not specified, then the function
%                        get host name by system('hostname') function.
%   - no_canlab          do not add canlab path
%   - no_cocoanlab       do not add cocoanlab path
%   - no_spm12           do not add spm12
%   - no_others          do not add other toolboxes (dicm2nii, BCT)
%   - ncpu               number of CPU cores to be used (default: 1).
%
%
% :Output:
% ::
%     rootdir            root directory based on hostname (not exactly
%                        'root' based on exact definition although...).
%     basedir            base working directory.
%     rscdir             resource directory that contains toolboxes.
%     gitdir             github directory that contains canlab and
%                        cocoanlab toolboxes.
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


[~, hostname] = system('hostname');
hostname = strtrim(hostname);
add_canlab = true;
add_cocoanlab = true;
add_spm12 = true;
add_others = true;
ncpu = 1;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'hostname'}
                hostname = varargin{i+1};
            case {'no_canlab'}
                add_canlab = false;
            case {'no_cocoanlab'}
                add_cocoanlab = false;
            case {'no_spm12'}
                add_spm12 = false;
            case {'no_others'}
                add_others = false;
            case {'ncpu'}
                ncpu = varargin{i+1};
        end
    end
end

setenv('DYLD_LIBRARY_PATH', [getenv('FREESURFER_HOME') '/lib/gcc/lib' ':/opt/X11/lib/flat_namespace']);
maxNumCompThreads(ncpu);

if regexp(hostname, ['\w*' 'hebenula.local' '\w*'])
    fprintf('Currently working on Habenula!\n');
    rootdir = '/Volumes/habenula';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
    gitdir = fullfile(rootdir, 'github');
elseif regexp(hostname, ['\w*' 'cocoanui-iMac-Pro4.local' '\w*'])
    fprintf('Currently working on alita!\n');
    rootdir = '/Volumes/alita';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
    gitdir = fullfile(rootdir, 'github');
elseif regexp(hostname, ['\w*' 'cnir' '\w*'])
    fprintf('Currently working on HPC!\n');
    rootdir = '/cocoanlab';
    basedir = fullfile(rootdir, 'habenula_sync2');
    rscdir = '/sas1/cocoanlab/Resources';
    gitdir = '/sas1/cocoanlab/Resources/github';
elseif regexp(hostname, ['\w*' 'JaeJoongui-MacBook-Pro.local' '\w*'])
    fprintf('Currently working on JJ macbook!\n');
    rootdir = '/Users/jaejoonglee';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
    gitdir = fullfile(rootdir, 'github');
else
    error('No matching hostname.');
end

if add_canlab; addpath(genpath(fullfile(gitdir, 'canlab'))); end
if add_cocoanlab; addpath(genpath(fullfile(gitdir, 'cocoanlab'))); end
if add_spm12; addpath(genpath(fullfile(rscdir, 'spm12'))); end
if add_others
    addpath(genpath(fullfile(rscdir, 'fmri_toolboxes/dicm2nii')));
    addpath(genpath(fullfile(rscdir, 'fmri_toolboxes/BCT')));
end


end