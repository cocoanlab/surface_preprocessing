function PREPROC = humanfmri_s12_smoothing(subject_code, study_imaging_dir, varargin)

% This function does smoothing of functional images using susan (FSL).
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s12_smoothing(subject_code, study_imaging_dir, varargin)
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
%   - fwhm               full-width half max for the smoothing kernel (mm).
%                        (default: 5 mm)
%
%
% :Output:
% ::
%     PREPROC.s_func_bold_files
%     PREPROC.s_func_reference_files
%     PREPROC.mean_s_func_bold_files
%     saves qc_images/sm_func_reference_masked.png 
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
fwhm = 5;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'fwhm'}
                fwhm = varargin{i+1};
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Smoothing', PREPROC.subject_code);
PREPROC.current_step = 's12';
PREPROC.current_step_letter = ['s' PREPROC.current_step_letter];

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        if PREPROC.ref_first_run
            ref_run = 1;
        else
            ref_run = i;
        end
        if regexp(PREPROC.current_step_letter, 'w')
            ref = PREPROC.w_func_reference_file{ref_run};
            ref_binmask = which('MNI152_T1_2mm_brain_mask.nii');
        else
            if regexp(PREPROC.current_step_letter, 'dc')
                ref = PREPROC.dc_func_reference_files{ref_run};
                ref_binmask = PREPROC.dc_func_reference_files_binarymask{ref_run};
            else
                if regexp(PREPROC.current_step_letter, 'r')
                    ref = PREPROC.func_reference_files{ref_run};
                    ref_binmask = PREPROC.func_reference_files_binarymask{ref_run};
                end
            end
        end
        
        if regexp(PREPROC.current_step_letter, 'w')
            input_dat = PREPROC.w_func_bold_files{i};
            mean_input_dat = PREPROC.mean_w_func_bold_files{i};
        else
            if regexp(PREPROC.current_step_letter, 'n')
                input_dat = PREPROC.n_func_bold_files{i};
                mean_input_dat = PREPROC.mean_n_func_bold_files{i};
            end
        end
        
        % Apply spatial smoothing
        fprintf('Applying spatial smoothing with FWHM %dmm kernel.\n', fwhm);
        PREPROC.s_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        [~, median_val] = system(['fslstats' ...
            ' ' mean_input_dat ...
            ' -k ' ref_mask ...
            ' -p 50']); % get median value
        median_val = str2num(median_val);
        brightness_threshold = median_val * 0.75;
        kernel_sigma = fwhm / sqrt(8*log(2));
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'susan' ...
            ' ' input_dat ...
            ' ' num2str(brightness_threshold) ...
            ' ' num2str(kernel_sigma) ...
            ' 3 1 1' ... % 3-dimension, use median, use mean image for edge detection
            ' ' mean_input_dat ...
            ' ' num2str(brightness_threshold) ...
            ' ' PREPROC.s_func_bold_files{i}]);
        system(['rm ' PREPROC.s_func_bold_files{i} '_usan_size.nii']);
        
        if ~ (PREPROC.ref_first_run && i ~= 1)
            
            fprintf('\n\nWorking on Reference image...\n\n');
            
            % Apply spatial smoothing
            fprintf('Applying spatial smoothing with FWHM %dmm kernel.\n', fwhm);
            PREPROC.s_func_reference_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_sm_reference.nii']);
            [~, median_val] = system(['fslstats' ...
                ' ' ref ...
                ' -k ' ref_mask ...
                ' -p 50']); % get median value
            median_val = str2num(median_val);
            brightness_threshold = median_val * 0.75;
            kernel_sigma = fwhm / sqrt(8*log(2));
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'susan' ...
                ' ' ref ...
                ' ' num2str(brightness_threshold) ...
                ' ' num2str(kernel_sigma) ...
                ' 3 1 0' ... % 3-dimension, use median, don't use other image for edge detection
                ' ' PREPROC.s_func_reference_files{i}]);
        end
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_s_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.s_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_s_func_bold_files{i}]);
        
    end
    
end

% Take snapshot of smoothed reference images
fprintf('Taking snapshot of smoothed reference images.\n');
canlab_preproc_show_montage(PREPROC.s_func_reference_files, fullfile(PREPROC.qcdir, 'sm_func_reference.png'));
drawnow;

% Take snapshot of smoothed mean functional images
fprintf('Taking snapshot of smoothed mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_s_func_bold_files, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold.png']));
drawnow;

close all;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end