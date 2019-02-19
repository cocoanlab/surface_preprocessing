function PREPROC = humanfmri_s2_surfrecon(subject_code, study_imaging_dir)

% This function does cortical surface reconstruction based on Freesurfer
% ('recon-all' function'), which is for further surfaced-based analysis.
% Please run the matlab through terminal with default environmental
% variable setting.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s2_surfrecon(subject_code, study_imaging_dir)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%
%
% :Output:
% ::
%     PREPROC.surfrecon_dir
%
%
% :Example:
% ::
%
%
% ..
%     Author and copyright information:
%
%     Copyright (C) Jan 2019  Jae-Joong Lee
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


fprintf('\n\n\n');
PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Surface Reconstruction', PREPROC.subject_code);
PREPROC.current_step = 's2';
PREPROC.current_step_letter = PREPROC.current_step_letter;

PREPROC.surfrecon_dir = fullfile(PREPROC.study_imaging_dir, 'surf');
if ~exist(PREPROC.surfrecon_dir, 'dir')
    mkdir(PREPROC.surfrecon_dir);
end

% Freesurfer work
system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
    ...
    'recon-all' ...
    ' -s ' PREPROC.subject_code ...
    ' -i ' PREPROC.anat_files{contains(PREPROC.anat_files, 'T1w')} ... % T1w
    ' -all']);

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end
