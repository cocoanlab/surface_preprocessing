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
anat_files = PREPROC.anat_files(~cellfun(@isempty, regexp(PREPROC.anat_files, 'T1w')));
n_anat = numel(anat_files);

if n_anat > 1
    % Longitudinal processing
    fprintf('Multiple anatomical images detected. Common anatomical template will be calculated.\n');
    [~, temp_surf_dir] = system('mktemp -d');
    temp_surf_dir = strtrim(temp_surf_dir);
    
    for i = 1:n_anat
        
        fprintf('Working on surface reconstruction: Timepoint %d ...\n', i);
        system(['export SUBJECTS_DIR=' temp_surf_dir ';' ...
            ...
            'recon-all' ...
            ' -s ' sprintf('TP%d', i) ...
            ' -i ' anat_files{i} ...
            ' -autorecon1' ...
            ' -gcareg' ...
            ' -canorm']);
        
    end
    
    fprintf('Working on surface reconstruction: Common template ...\n');
    system(['export SUBJECTS_DIR=' temp_surf_dir ';' ...
        ...
        'recon-all' ...
        ' -base ' PREPROC.subject_code ...
        reshape(strcat(repmat(' -tp TP', n_anat, 1), num2str([1:n_anat]'))', 1, []) ...
        ' -all']);
    
    system(['mv ' fullfile(temp_surf_dir, PREPROC.subject_code) ' ' PREPROC.surfrecon_dir]);
    system(['rm ' temp_surf_dir]);
    
else

    fprintf('Working on surface reconstruction ...\n');
    system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
        ...
        'recon-all' ...
        ' -s ' PREPROC.subject_code ...
        ' -i ' anat_files{1} ...
        ' -all']);

end

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end
