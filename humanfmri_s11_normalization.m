function PREPROC = humanfmri_s11_normalization(subject_code, study_imaging_dir, varargin)

% This function does normalization of functional images to MNI space.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s11_normalization(subject_code, study_imaging_dir, varargin)
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
%   - run_num            runs to include. ex) [1 2 4 5], ...
%   - n_thread            the number of threads for ANTs (Default: 1)
%
%
% :Output:
% ::
%     PREPROC.w_func_bold_files
%     PREPROC.w_func_reference_file
%     PREPROC.mean_w_func_bold_files
%     saves qc_images/mean_[prefix]_func_bold_masked.png
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
do_select_run = false;
run_num = NaN;
n_thread = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'n_thread'}
                n_thread = varargin{i+1};
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Normalization', PREPROC.subject_code);
PREPROC.current_step = 's11';
PREPROC.current_step_letter = ['w' PREPROC.current_step_letter];

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        % Warp functional BOLD images to MNI space
        fprintf('Warping functional BOLD images to MNI space...\n');
        PREPROC.w_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        
        if isfield(PREPROC, 'surfrecon_dir')
            [~, volinfo] = system(['fslnvols ' PREPROC.n_func_bold_files{i}]);
            num_vols = str2num(volinfo);
            [~, temp_prewarp_dir] = system('mktemp -d');
            temp_prewarp_dir = strtrim(temp_prewarp_dir);
            [~, temp_postwarp_dir] = system('mktemp -d');
            temp_postwarp_dir = strtrim(temp_postwarp_dir);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'fslsplit' ...
                ' ' PREPROC.n_func_bold_files{i} ...
                ' ' fullfile(temp_prewarp_dir, 'prevol')]);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'for i in {0..' num2str(PREPROC.func_nvol(i)-1) '}; do ' ...
                'applywarp' ...
                ' --in=' fullfile(temp_prewarp_dir, 'prevol') '$(printf "%.4d" $i).nii' ...
                ' --ref=' PREPROC.anat_reference_file_warped ...
                ' --warp=' PREPROC.anat_reference_file_warpfield ... % T1 to MNI
                ' --premat=' PREPROC.transform_EPI_to_T1 ... % EPI to T1
                ' --interp=spline' ...
                ' --out=' fullfile(temp_postwarp_dir, 'postvol') '$(printf "%.4d" $i).nii' ...
                '; done']);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'cd ' temp_postwarp_dir ';' ...
                'fslmerge' ...
                ' -tr' ...
                ' ' PREPROC.w_func_bold_files{i} ...
                ' ' 'postvol*.nii' ...
                ' ' num2str(PREPROC.func_TR(i) / 1000)]); % sec
            system(['rm -r ' temp_prewarp_dir]);
            system(['rm -r ' temp_postwarp_dir]);
        else
            system(['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=' num2str(n_thread) ';' ...
                'antsApplyTransforms' ...
                ' --verbose 1' ...
                ' --dimensionality 3' ...
                ' --input-image-type 3' ...
                ' --float 1' ...
                ' --input ' PREPROC.n_func_bold_files{i} ...
                ' --interpolation BSpline' ...
                ' --output ' PREPROC.w_func_bold_files{i} ...
                ' --reference-image ' PREPROC.anat_reference_file_warped ...
                ' --transform ' PREPROC.anat_reference_file_warpfield ...
                ' --transform ' PREPROC.anat_reference_file_warpmat ...
                ' --transform ' PREPROC.transform_EPI_to_T1_ITK]);
        end
        
        if i == 1
            fprintf('\n\nWorking on Reference image...\n\n');
            
            % Warp functional reference image to MNI space
            fprintf('Warping functional reference image to MNI space...\n');
            PREPROC.w_func_reference_file = fullfile(PREPROC.preproc_func_dir, 'normalized_to_MNI_func_reference.nii');
%             PREPROC.w_func_reference_file_masked = fullfile(PREPROC.preproc_func_dir, 'normalized_to_MNI_func_reference_masked.nii');
            if contains(PREPROC.current_step_letter, 'dc')
                func_ref = PREPROC.dc_func_reference_file_masked;
            else
                if contains(PREPROC.current_step_letter, 'r')
                    func_ref = PREPROC.func_reference_file_masked;
                end
            end
            if isfield(PREPROC, 'surfrecon_dir')
                system(['export FSLOUTPUTTYPE=NIFTI;' ...
                    ...
                    'applywarp' ...
                    ' --in=' func_ref ...
                    ' --ref=' PREPROC.anat_reference_file_warped ...
                    ' --warp=' PREPROC.anat_reference_file_warpfield ...
                    ' --premat=' PREPROC.transform_EPI_to_T1 ...
                    ' --interp=spline' ...
                    ' --out=' PREPROC.w_func_reference_file]);
%                 system(['export FSLOUTPUTTYPE=NIFTI;' ...
%                     ...
%                     'fslmaths' ...
%                     ' ' PREPROC.w_func_reference_file ...
%                     ' -mul ' PREPROC.anat_reference_file_warped_binarymask ...
%                     ' ' PREPROC.w_func_reference_file_masked]);
            else
                system(['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=' num2str(n_thread) ';' ...
                    'antsApplyTransforms' ...
                    ' --verbose 1' ...
                    ' --dimensionality 3' ...
                    ' --input-image-type 0' ...
                    ' --float 1' ...
                    ' --input ' func_ref ...
                    ' --interpolation BSpline' ...
                    ' --output ' PREPROC.w_func_reference_file ...
                    ' --reference-image ' PREPROC.anat_reference_file_warped ...
                    ' --transform ' PREPROC.anat_reference_file_warpfield ...
                    ' --transform ' PREPROC.anat_reference_file_warpmat ...
                    ' --transform ' PREPROC.transform_EPI_to_T1_ITK]);
            end
        end
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_w_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.w_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_w_func_bold_files{i}]);
        
    end
    
end

% Take snapshot of warped mean functional images
fprintf('Taking snapshot of warped mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_w_func_bold_files, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold.png']));
drawnow;

close all;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end