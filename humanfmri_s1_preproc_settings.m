function PREPROC = humanfmri_s1_preproc_settings(subject_code, study_imaging_dir, varargin)

% This function sets up initial environment for data preprocessing such as
% creating directories, choosing reference image, and save mean functional
% image for quality-check.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s1_preproc_settings(subject_code, study_imaging_dir, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%
%
% :Optional Input:
% ::
%   - no_update          do not update PREPROC.mat file
%   - no_fmap            do not make fieldmap folder
%   - run_num            runs to include. ex) [1 2 4 5], ...
%   - ref_first_run      use reference image of first run for whole
%                        preprocessing. (default: false)
%
%
% :Output:
% ::
%     PREPROC.current_step
%     PREPROC.current_step_letter
%     PREPROC.ref_first_run
%     PREPROC.preproc_outputdir
%     PREPROC.preproc_anat_dir
%     PREPROC.preproc_func_dir
%     PREPROC.preproc_mean_func_dir
%     PREPROC.preproc_fmap_dir (if selected)
%     PREPROC.qcdir
%     PREPROC.func_reference_files
%     PREPROC.func_reference_files_masked
%     PREPROC.func_reference_files_binarymask
%     PREPROC.setting_reference_image
%     PREPROC.mean_func_bold_files
%     PREPROC.mean_func_bold_files_masked
%     saves qc_images/func_reference_masked.png 
%     saves qc_images/mean_func_bold_masked.png
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


fprintf('\n\n\n');
do_update = true;
make_fmap = true;
do_select_run = false;
run_num = NaN;
ref_first_run = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'no_update'}
                do_update = false;
            case {'no_fmap'}
                make_fmap = false;
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'ref_first_run'}
                ref_first_run = true;
        end
    end
end


subject_dir = fullfile(study_imaging_dir, 'preprocessed', subject_code);
if do_update && exist(fullfile(subject_dir, 'PREPROC.mat')) == 2 % do update & mat file exist
    PREPROC = save_load_PREPROC(subject_dir, 'load'); % load PREPROC
else
    preproc_mat_dir = fullfile(study_imaging_dir, 'raw', subject_code);
    PREPROC = save_load_PREPROC(preproc_mat_dir, 'load'); % load PREPROC
end

print_header('Setting up initial environment for preprocessing', PREPROC.subject_code);
PREPROC.current_step = 's1';
PREPROC.current_step_letter = '';

PREPROC.ref_first_run = ref_first_run;

fprintf('Making directories for preprocessing.\n');
PREPROC.preproc_outputdir = subject_dir;
system(['mkdir -p ' PREPROC.preproc_outputdir]);
PREPROC.preproc_anat_dir = fullfile(PREPROC.preproc_outputdir, 'anat');
system(['mkdir -p ' PREPROC.preproc_anat_dir]);
PREPROC.preproc_func_dir = fullfile(PREPROC.preproc_outputdir, 'func');
system(['mkdir -p ' PREPROC.preproc_func_dir]);
PREPROC.preproc_mean_func_dir = fullfile(PREPROC.preproc_outputdir, 'mean_func');
system(['mkdir -p ' PREPROC.preproc_mean_func_dir]);
if make_fmap
    PREPROC.preproc_fmap_dir = fullfile(PREPROC.preproc_outputdir, 'fmap');
    system(['mkdir -p ' PREPROC.preproc_fmap_dir]);
end
PREPROC.qcdir = fullfile(PREPROC.preproc_outputdir, 'qc_images');
system(['mkdir -p ' PREPROC.qcdir]);

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        if ~ (PREPROC.ref_first_run && i ~= 1)
            
            % Select reference image
            PREPROC.func_reference_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_reference.nii']);
            if exist(PREPROC.func_sbref_files{i}) % if SBRef exists, it is used as a reference image
                fprintf('Use SBRef image as reference.\n');
                PREPROC.setting_reference_image{i, 1} = 'sbref';
                system(['cp' ...
                    ' ' PREPROC.func_sbref_files{i} ...
                    ' ' PREPROC.func_reference_files{i}]); % Use SBRef image
            else  % if SBRef does not exist, first volume of BOLD image is used as a reference image
                fprintf('Use first volume of BOLD image as reference.\n');
                PREPROC.setting_reference_image{i, 1} = 'first_vol';
                system(['export FSLOUTPUTTYPE=NIFTI;' ...
                    ...
                    'fslroi' ...
                    ' ' PREPROC.func_bold_files{i} ...
                    ' ' PREPROC.func_reference_files{i} ...
                    ' 0 1']); % Use first volume
            end
            
            % Extract brain mask of reference image
            fprintf('Extracting brain of reference image and saving implicit mask...\n');
            PREPROC.func_reference_files_masked{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_reference_masked.nii']);
            PREPROC.func_reference_files_binarymask{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_reference_masked_mask.nii']);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'bet' ...
                ' ' PREPROC.func_reference_files{i} ...
                ' ' PREPROC.func_reference_files_masked{i} ...
                ' -f 0.3' ... % ICA-AROMA recommendation
                ' -n' ...
                ' -m' ...
                ' -R']);
            
        end
        
        if PREPROC.ref_first_run
            ref_run = 1;
        else
            ref_run = i;
        end
        ref = PREPROC.func_reference_files{ref_run};
        ref_masked = PREPROC.func_reference_files_masked{ref_run};
        ref_binmask = PREPROC.func_reference_files_binarymask{ref_run};
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' b '.nii']);
        PREPROC.mean_func_bold_files_masked{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' b '_masked.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_func_bold_files{i}]);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.mean_func_bold_files{i} ...
            ' -mul ' ref_binmask ...
            ' ' PREPROC.mean_func_bold_files_masked{i}]);
        
    end
    
end

% Take snapshot of reference images
fprintf('Taking snapshot of reference images.\n');
canlab_preproc_show_montage(PREPROC.func_reference_files_masked, fullfile(PREPROC.qcdir, 'func_reference_masked.png'));
drawnow;

% Take snapshot of mean functional images
fprintf('Taking snapshot of mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_func_bold_files_masked, fullfile(PREPROC.qcdir, 'mean_func_bold_masked.png'));
drawnow;

close all;

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % load PREPROC
fprintf('\n\n\n');

end