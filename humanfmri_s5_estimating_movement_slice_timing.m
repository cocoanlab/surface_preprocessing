function PREPROC = humanfmri_s5_estimating_movement_slice_timing(subject_code, study_imaging_dir, varargin)

% This function does slice time correction for functional data using 
% 3dTshift (AFNI) function, with extracting slice timing information from
% "MosaicRefAcqTimes" provided by dicm2nii (DICOM header file).
% Since slice timing correction hampers acquisition of exact realignment
% parameters (movement), estimating realignment parameters is preceded
% before slice timing correction using 3dvolreg (AFNI) function without any
% output images.
% Note: Skip this step for multi-band dataset. 
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s5_estimating_movement_slice_timing(subject_code, study_imaging_dir, varargin)
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
%   - custom_timing      specification of slice timing. order(integer form)
%                        or timing (in milliseconds, float form).
%                        (e.g., order : {[1 3 5 7 2 4 6 8]} [for all run]
%                               timing : {[0.0000 252.5000 62.5000
%                                        315.0000 125.0000 377.5000
%                                        190.0000]} [for all run])
%                        to specify different slice timing for each run,
%                        just add matrix in cell variable.
%                        (e.g., order : {[1 3 5 2 4], [2 4 1 3 5], ...})
%                        if not specified, it is obtained by reading dicom
%                        header file.
%
%
% :Output:
% ::
%     PREPROC.mvmt_param_files
%     PREPROC.mvmt_matrix_files
%     PREPROC.slice_timing_files
%     PREPROC.a_func_bold_files
%     PREPROC.mean_a_func_bold_files
%     PREPROC.mean_a_func_bold_files_masked
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
custom_timing = {[]};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'custom_timing'}
                custom_timing = varargin{i+1};
                if ~iscell(custom_timing)
                    custom_timing = {custom_timing};
                end
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Estimating movement & Slice timing correction', PREPROC.subject_code);
PREPROC.current_step = 's5';
PREPROC.current_step_letter = ['a' PREPROC.current_step_letter];

if numel(custom_timing) == 1
    custom_timing = repmat(custom_timing, numel(PREPROC.func_bold_files), 1);
end
custom_timing = custom_timing(:);

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        % Estimate movement parameter
        fprintf('Estimating movement parameter before slice timing correction...\n');
        PREPROC.mvmt_param_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['rp_' b '.1D']);
        PREPROC.mvmt_matrix_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['rp_mat_' b '.1D']);
        system(['3dvolreg' ...
            ' -base ' PREPROC.func_reference_file ...
            ' -verbose' ...
            ' -Fourier' ...
            ' -twopass' ...
            ' -zpad 4' ...
            ' -1Dfile ' PREPROC.mvmt_param_files{i} ... % Roll/Pitch/Yaw (deg, counter-clockwise) & dS/dL/dP (mm)
            ' -1Dmatrix_save ' PREPROC.mvmt_matrix_files{i} ...
            ' -float' ...
            ' -prefix NULL' ...
            ' ' PREPROC.func_bold_files{i}]);
        
        % Slice timing correction
        fprintf('Slice timing correction...\n');
        if isempty(custom_timing{i})
            load(PREPROC.func_bold_dicomheader_files{1});
            custom_timing{i} = h.MosaicRefAcqTimes;
        else
            if numel(custom_timing{i}) ~= PREPROC.func_nslices(i)
                error('Number of slice timing does not match with data. Please check the exact slice timing.');
            end
            if min(custom_timing{i}) == 1 && all(diff(unique(custom_timing{i})) == 1) % Slice order was supplied
                custom_timing{i} = (custom_timing{i} - 1) * PREPROC.func_TR(i) / max(custom_timing{i});
            elseif min(custom_timing{i}) == 0 && max(custom_timing{i}) < PREPROC.func_TR(i) % Slice timing was supplied
                custom_timing{i} = custom_timing{i};
            else
                error('Wrong custom slice timing. Please check the exact slice timing.')
            end
        end
        custom_timing{i} = custom_timing{i}(:); % Vectorization!
        
        fprintf('\nSlice timing information:\n\n');
        disp(custom_timing{i});
        PREPROC.slice_timing_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['slice_timing_' b '.txt']);
        fid = fopen(PREPROC.slice_timing_files{i}, 'wt');
        fprintf(fid, '%.4f\t', custom_timing{i});
        fclose(fid);

        PREPROC.a_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        system(['3dTshift' ...
            ' -verbose' ...
            ' -Fourier' ...
            ' -TR ' num2str(PREPROC.func_TR(i)) 'ms' ...
            ' -tpattern @' PREPROC.slice_timing_files{i} ...
            ' -prefix ' PREPROC.a_func_bold_files{i} ...
            ' ' PREPROC.func_bold_files{i}]);
        %             ' -tzero ' '0' ... <- Reference slice as first! Deafult is mean slice
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_a_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        PREPROC.mean_a_func_bold_files_masked{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '_masked.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.a_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_a_func_bold_files{i}]);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.mean_a_func_bold_files{i} ...
            ' -mul ' PREPROC.func_reference_file_binarymask ...
            ' ' PREPROC.mean_a_func_bold_files_masked{i}]);
        
    end
    
end

% Take snapshot of slice timing-corrected mean functional images
fprintf('Taking snapshot of slice timing-corrected mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_a_func_bold_files_masked, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold_masked.png']));
drawnow;

close all;

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end


