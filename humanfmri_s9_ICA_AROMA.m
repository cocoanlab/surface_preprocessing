function PREPROC = humanfmri_s9_ICA_AROMA(subject_code, study_imaging_dir, ica_aroma_basedir, varargin)

% This function de-noises movement related signals using ICA-AROMA.
% Before the process, smooothing is temporaily applied to get robust ICs.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s9_ICA_AROMA(subject_code, study_imaging_dir, ica_aroma_basedir, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%   - ica_aroma_basedir  the ICA-AROMA toolbox directory.
%
%
% :Optional Input:
% ::
%   - run_num            runs to include. ex) [1 2 4 5], ...
%   - fwhm               full-width half max for the smoothing kernel (mm).
%                        (default: 5 mm)
%   - n_dim              sets dimensionality for melodic (number of ICs).
%                        (e.g., {'200'} [for all run]
%                               or {'200', '100', '200'} [for each run])
%                        if not specified, it is automatically estimated by
%                        melodic algorithm.
%   - n_thread           the number of threads for ANTs (Default: 1)
%   - filt_on_sm         Filter noise-ICs on smoothed data (Default: false)
%
%
% :Output:
% ::
%     PREPROC.ica_aroma_dir
%     PREPROC.i_func_bold_files
%     PREPROC.mean_i_func_bold_files
%     PREPROC.mean_i_func_bold_files_masked
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
n_dim = {0};
n_thread = 1;
filt_on_sm = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'fwhm'}
                fwhm = varargin{i+1};
            case {'n_dim'}
                n_dim = varargin{i+1};
                if ~iscell(n_dim)
                    n_dim = {n_dim};
                end
            case {'n_thread'}
                n_thread = varargin{i+1};
            case {'filt_on_sm'}
                filt_on_sm = true;
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('ICA-AROMA', PREPROC.subject_code);
PREPROC.current_step = 's9';
if filt_on_sm
    PREPROC.current_step_letter = ['is' PREPROC.current_step_letter];
else
    PREPROC.current_step_letter = ['i' PREPROC.current_step_letter];
end

if numel(n_dim) == 1
    n_dim = repmat(n_dim, numel(PREPROC.func_bold_files), 1);
end
n_dim = n_dim(:);

if regexp(PREPROC.current_step_letter, 'w')
    func_ref_mask = '$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz';
else
    if regexp(PREPROC.current_step_letter, 'dc')
        func_ref_mask = PREPROC.dc_func_reference_file_binarymask;
    else
        if regexp(PREPROC.current_step_letter, 'r')
            func_ref_mask = PREPROC.func_reference_file_binarymask;
        end
    end
end

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        if regexp(PREPROC.current_step_letter, 's')
            
            input_dat = PREPROC.s_func_bold_files{i};
            mean_input_dat = PREPROC.mean_s_func_bold_files{i};
            input_denoising = input_dat;
            
        else
            
            if regexp(PREPROC.current_step_letter, 'w')
                unsmoothed_input_dat = PREPROC.w_func_bold_files{i};
                unsmoothed_mean_input_dat = PREPROC.mean_w_func_bold_files{i};
            else
                if regexp(PREPROC.current_step_letter, 'dc')
                    unsmoothed_input_dat = PREPROC.dc_func_bold_files{i};
                    unsmoothed_mean_input_dat = PREPROC.mean_dc_func_bold_files{i};
                else
                    if regexp(PREPROC.current_step_letter, 'r')
                        unsmoothed_input_dat = PREPROC.r_func_bold_files{i};
                        unsmoothed_mean_input_dat = PREPROC.mean_r_func_bold_files{i};
                    end
                end
            end
            
            % Apply spatial smoothing for extracting good ICs
            fprintf('Applying spatial smoothing for extracting good ICs...\n');
            [~, temp_sm_dir] = system('mktemp -d');
            temp_sm_dir = strtrim(temp_sm_dir);
            input_dat = fullfile(temp_sm_dir, ['temp_sm_' b '.nii']);
            [~, median_val] = system(['fslstats' ...
                ' ' unsmoothed_mean_input_dat ...
                ' -k ' func_ref_mask ...
                ' -p 50']); % get median value
            median_val = str2num(median_val);
            brightness_threshold = median_val * 0.75;
            kernel_sigma = fwhm / sqrt(8*log(2));
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'susan' ...
                ' ' unsmoothed_input_dat ...
                ' ' num2str(brightness_threshold) ...
                ' ' num2str(kernel_sigma) ...
                ' 3 1 1' ... % 3-dimension, use median, use mean image for edge detection
                ' ' unsmoothed_mean_input_dat ...
                ' ' num2str(brightness_threshold) ...
                ' ' input_dat]);
            
            if filt_on_sm
                input_denoising = input_dat;
            else
                input_denoising = unsmoothed_input_dat;
            end
            
        end
        
        % ICA-AROMA
        fprintf('Running ICA-AROMA...');
        [~, functype] = fileparts(PREPROC.dicom_func_bold_dir{i});
        PREPROC.ica_aroma_dir{i, 1} = fullfile(PREPROC.preproc_func_dir, ['ica_aroma_' functype]);
        
        if isfield(PREPROC, 'surfrecon_dir')
            system(['export FSLOUTPUTTYPE=NIFTI_GZ;' ...
                ...
                'python2.7 ' fullfile(ica_aroma_basedir, 'ICA_AROMA.py') ...
                ' -i ' input_dat ...
                ' -a ' PREPROC.transform_EPI_to_T1 ...
                ' -w ' PREPROC.anat_reference_file_warpfield ...
                ' -mc ' PREPROC.mvmt_param_files_demeaned{i} ...
                ' -m ' func_ref_mask ...
                ' -tr ' num2str(PREPROC.func_TR(i) / 1000) ... % sec
                ' -dim ' num2str(n_dim{i}) ...
                ' -den ' 'no' ...
                ' -o ' PREPROC.ica_aroma_dir{i}]);
        else
            system(['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=' num2str(n_thread) ';' ...
                'export FSLOUTPUTTYPE=NIFTI_GZ;' ...
                ...
                'python2.7 ' fullfile(ica_aroma_basedir, 'ICA_AROMA.py') ...
                ' -i ' input_dat ...
                ' -ITKcoreg ' PREPROC.transform_EPI_to_T1_ITK ...
                ' -ITKaffine ' PREPROC.anat_reference_file_warpmat ...
                ' -ITKwarp ' PREPROC.anat_reference_file_warpfield ...
                ' -mc ' PREPROC.mvmt_param_files_demeaned{i} ...
                ' -m ' func_ref_mask ...
                ' -tr ' num2str(PREPROC.func_TR(i) / 1000) ... % sec
                ' -dim ' num2str(n_dim{i}) ...
                ' -den ' 'no' ...
                ' -o ' PREPROC.ica_aroma_dir{i}]);
        end
        motion_ICs = importdata(fullfile(PREPROC.ica_aroma_dir{i}, 'classified_motion_ICs.txt'));
        motion_ICs_str = regexprep(num2str(motion_ICs, '%d,'), ',\s+', ',');
        motion_ICs_str = ['"' motion_ICs_str(1:end-1) '"'];
        PREPROC.i_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fsl_regfilt' ...
            ' -v' ...
            ' --in=' input_denoising ...
            ' --design=' fullfile(PREPROC.ica_aroma_dir{i}, 'melodic.ica', 'melodic_mix') ...
            ' --mask=' func_ref_mask ... % Output will be masked!!!
            ' --filter=' motion_ICs_str ...
            ' --out=' PREPROC.i_func_bold_files{i}]);
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_i_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.i_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_i_func_bold_files{i}]);
        
        if exist('temp_sm_dir', 'var')
            system(['rm -r ' temp_sm_dir]);
        end
        
    end
    
end

% Take snapshot of ICA-AROMA-denoised mean functional images
fprintf('Taking snapshot of ICA-AROMA-denoised mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_i_func_bold_files, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold.png']));
drawnow;

close all;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end