function [rootdir, basedir, rscdir] = set_path_env(varargin)

% The function gets useful directory names and sets basic environment
% (path to toolboxes, environment variables, maximum number of CPU cores)
% depending on current hostname.
%
%
% :Usage:
% ::
%     [rootdir, basedir, rscdir] = set_path_env(varargin)
%
%
% :Optional Input:
% ::
%   - hostname           custom hostname.
%                        if this option is not specified, then the function
%                        get host name by system('hostname') function.
%   - ncpu               number of CPU cores to be used (default: 1).
%
%
% :Output:
% ::
%     rootdir            root directory based on hostname (not exactly
%                        'root' based on exact definition although...).
%     basedir            base working directory.
%     rscdir             resource directory that contains toolboxes.
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
ncpu = 1;
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'hostname'}
                hostname = varargin{i+1};
            case {'ncpu'}
                ncpu = varargin{i+1};
        end
    end
end

if regexp(hostname, ['\w*' 'hebenula.local' '\w*'])
    fprintf('Currently working on Habenula!\n');
    rootdir = '/Volumes/habenula';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
elseif regexp(hostname, ['\w*' 'cocoanui-iMac-Pro4.local' '\w*'])
    fprintf('Currently working on alita!\n');
    rootdir = '/Volumes/alita';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
elseif regexp(hostname, ['\w*' 'cnir' '\w*'])
    fprintf('Currently working on HPC!\n');
    rootdir = '/cocoanlab';
    basedir = fullfile(rootdir, 'habenula_sync2');
    rscdir = fullfile(rootdir, 'Resources');
elseif regexp(hostname, ['\w*' 'JaeJoongui-MacBook-Pro.local' '\w*'])
    fprintf('Currently working on JJ macbook!\n');
    rootdir = '/Users/jaejoonglee';
    basedir = fullfile(rootdir, 'hbmnas');
    rscdir = fullfile(rootdir, 'Resource');
else
    error('No matching hostname.');
end

addpath(genpath(fullfile(rscdir, 'matlab_toolboxes/canlab')));
addpath(genpath(fullfile(rscdir, 'matlab_toolboxes/cocoanlab')));
addpath(genpath(fullfile(rscdir, 'spm12')));

rmpath(genpath(fullfile(rscdir, 'matlab_toolboxes/canlab/CanlabPrivate/dicm2nii')));
rmpath(genpath(fullfile(rscdir, 'matlab_toolboxes/cocoanlab/humanfmri_preproc_bids/external/dicm2nii')));
addpath(genpath(fullfile(rscdir, 'matlab_toolboxes/dicm2nii')));

setenv('DYLD_LIBRARY_PATH', [getenv('FREESURFER_HOME') '/lib/gcc/lib' ':/opt/X11/lib/flat_namespace']);
maxNumCompThreads(ncpu);

end