function PREPROC = humanfmri_s7_distortion_correction(subject_code, study_imaging_dir, varargin)

% This function applies the distortion correction using TOPUP and
% APPLYTOPUP function (FSL). Encoding direction will be extracted from
% "UnwarpDirection" provided by dicm2nii (DICOM header file).
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s7_distortion_correction(subject_code, study_imaging_dir)
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
%
%   - run_num            runs to include. ex) [1 2 4 5], ...
%   - custom_acqparam    specification of acquisition parameters for
%                        distortion correction.
%                        (e.g.,    encoding dir    readout time
%                              R->L   P->A   I->S     (sec)
%                               [0     -1      0      0.0627
%                                0     -1      0      0.0627
%                                0      1      0      0.0627
%                                0      1      0      0.0627])
%                        if not specified, these parameters are obtained by
%                        reading fieldmap image and dicom header file.
%   - epi_enc_dir        Phase encoding direction of functional data.
%                        (e.g., 'ap'
%                               -> direction is A->P for every run,
%                               {'ap', 'pa', 'lr'}
%                               -> direction is A->P for run 1,
%                               -> direction is P->A for run 2,
%                               -> direction is L->R for run 3)
%                        if not specified, it is obtained by reading dicom
%                        header file.
%
%
% :Output:
% ::   
%     PREPROC.topup.fmap_combined_file
%                  .dc_param_file
%                  .output_path
%                  .fieldmap_file
%                  .unwarped_file
%                  .config_file
%     PREPROC.dc_func_bold_files
%     PREPROC.dc_func_reference_files
%     PREPROC.dc_func_reference_files_masked
%     PREPROC.dc_func_reference_files_binarymask
%     PREPROC.mean_dc_func_bold_files
%     PREPROC.mean_dc_func_bold_files_masked
%     saves qc_images/dc_func_reference_masked.png 
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
dc_param_dat = [];
epi_enc_dir = {''};

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'custom_acqparam'}
                dc_param_dat = varargin{i+1};
            case {'epi_enc_dir'}
                epi_enc_dir = varargin{i+1};
                if ~iscell(epi_enc_dir)
                    epi_enc_dir = {epi_enc_dir};
                end
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Distortion correction', PREPROC.subject_code);
PREPROC.current_step = 's7';
PREPROC.current_step_letter = ['dc' PREPROC.current_step_letter];

if numel(epi_enc_dir) == 1
    epi_enc_dir = repmat(epi_enc_dir, numel(PREPROC.func_bold_files), 1);
end
epi_enc_dir = epi_enc_dir(:);

fmap_enc_dictionary = ... % LAS system! (radiological convention)
    {'x', 'rl', [1 0 0]; ...
    '-x', 'lr', [-1 0 0]; ...
    'y', 'pa', [0 1 0]; ...
    '-y', 'ap', [0 -1 0]; ...
    'z', 'is', [0 0 1]; ...
    '-z', 'si', [0 0 -1]};

if ~isfield(PREPROC, 'topup')
   
    % Combine fieldmap images
    fprintf('Combining fieldmap images.\n');
    PREPROC.topup.fmap_combined_file = fullfile(PREPROC.preproc_fmap_dir, [PREPROC.subject_code '_dir-combined_epi.nii']);
    fmap_filelist = strcat(PREPROC.fmap_files, {' '});
    fmap_filelist = cat(2, fmap_filelist{:});
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'fslmerge' ...
        ' -t' ...
        ' ' PREPROC.topup.fmap_combined_file ...
        ' ' fmap_filelist]);

    % Extract acqusition parameters for TOPUP
    fprintf('Extracting acqusition parameters for TOPUP.\n');

    if isempty(dc_param_dat)

        for i = 1:numel(PREPROC.fmap_files)

            load(PREPROC.fmap_dicomheader_files{i});

            [regexp_start, regexp_end] = regexp(PREPROC.fmap_files{i}, 'dir-\w\w_epi');
            fmap_enc_dir_label = PREPROC.fmap_files{i}(regexp_start+length('dir-'):regexp_end-length('_epi'));
            if strcmp(fmap_enc_dictionary(:,1), h.UnwarpDirection) == strcmp(fmap_enc_dictionary(:,2), fmap_enc_dir_label)
                fmap_enc_dir = fmap_enc_dictionary{strcmp(fmap_enc_dictionary(:,1), h.UnwarpDirection), 3};
            else
                error('Unidentifiable phase encoding direction of fieldmap images.')
            end

            readout_time = h.ReadoutSeconds;
            distort_num = h.NumberOfTemporalPositions;
            dc_param_dat = [dc_param_dat; repmat([fmap_enc_dir readout_time], distort_num, 1)];

        end

    end

    PREPROC.topup.dc_param_file = fullfile(PREPROC.preproc_fmap_dir, ['dc_param_', PREPROC.subject_code '_combined_epi.txt']);
    fid = fopen(PREPROC.topup.dc_param_file, 'wt');
    fprintf(fid, [repmat('%.4f\t', 1, size(dc_param_dat, 2)), '\n'], dc_param_dat');
    fclose(fid);
    fprintf('\n\nAcqusition parameters for distortion correction:\n\n');
    disp(array2table(dc_param_dat, 'VariableNames', {'R_L', 'P_A', 'I_S', 'ReadoutTime'}));

    % Run TOPUP
    fprintf('Running TOPUP to estimate distortion...\n');
    PREPROC.topup.output_path = fullfile(PREPROC.preproc_fmap_dir, ['topup_output_', PREPROC.subject_code]);
    PREPROC.topup.fieldmap_file = fullfile(PREPROC.preproc_fmap_dir, ['topup_fieldmap_', PREPROC.subject_code '.nii']);
    PREPROC.topup.unwarped_file = fullfile(PREPROC.preproc_fmap_dir, ['topup_unwarped_', PREPROC.subject_code '.nii']);
    PREPROC.topup.config_file = [getenv('FSLDIR') '/etc/flirtsch/b02b0.cnf'];
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'topup' ...
        ' --verbose' ...
        ' --imain=' PREPROC.topup.fmap_combined_file ...
        ' --datain=' PREPROC.topup.dc_param_file ...
        ' --config=' PREPROC.topup.config_file ...
        ' --out=' PREPROC.topup.output_path ...
        ' --fout=' PREPROC.topup.fieldmap_file ...
        ' --iout=' PREPROC.topup.unwarped_file]);

    % Take snapshot of fieldmap images before/after distortion correction
    fprintf('Take snapshot of fieldmap images before/after TOPUP.\n');
    j = 0;
    for i = 1:numel(PREPROC.fmap_files)
        [~, volinfo] = system(['echo $(fslval ' PREPROC.fmap_files{i} ' dim4)']);
        num_vols = str2num(volinfo);
        fmap_idx = (j+1):(j+num_vols);
        [~, fmap_name] = fileparts(PREPROC.fmap_files{i});
        topup_unwarped_png = fullfile(PREPROC.qcdir, ['topup_unwarped_' fmap_name '.png']);
        topup_before_list = cellstr(strcat(PREPROC.topup.fmap_combined_file, ',', num2str(fmap_idx')));
        topup_after_list = cellstr(strcat(PREPROC.topup.unwarped_file, ',', num2str(fmap_idx')));
        canlab_preproc_show_montage([topup_before_list; topup_after_list], topup_unwarped_png);
        drawnow;
        close all;
        j = max(fmap_idx);
    end
    
else
    
    fprintf('Existing TOPUP files found.\n');
    dc_param_dat = importdata(PREPROC.topup.dc_param_file);
    fprintf('\n\nAcqusition parameters for distortion correction:\n\n');
    disp(array2table(dc_param_dat, 'VariableNames', {'R_L', 'P_A', 'I_S', 'ReadoutTime'}));
    
end

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        % Run APPLYTOPUP
        fprintf('Correcting distortion by APPLYTOPUP...\n');
        PREPROC.dc_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        if isempty(epi_enc_dir{i})
            load(PREPROC.func_bold_dicomheader_files{i});
            applytopup_idx = find(ismember(dc_param_dat(:, 1:3), fmap_enc_dictionary{strcmp(fmap_enc_dictionary(:,1), h.UnwarpDirection), 3}, 'rows'), 1); % First index only
        else
            applytopup_idx = find(ismember(dc_param_dat(:, 1:3), fmap_enc_dictionary{strcmp(fmap_enc_dictionary(:,2), epi_enc_dir{i}), 3}, 'rows'), 1); % First index only
        end
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'applytopup' ...
            ' --verbose' ...
            ' --imain=' PREPROC.r_func_bold_files{i} ...
            ' --inindex=' num2str(applytopup_idx) ...
            ' --topup=' PREPROC.topup.output_path ...
            ' --datain=' PREPROC.topup.dc_param_file ...
            ' --method=jac' ...
            ' --interp=spline' ...
            ' --out=' PREPROC.dc_func_bold_files{i}]);
        
        % Removing spline interpolation neg values by absolute
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.dc_func_bold_files{i} ...
            ' -abs' ...
            ' ' PREPROC.dc_func_bold_files{i}]);
        
        if ~ (PREPROC.ref_first_run && i ~= 1)
            
            fprintf('\n\nWorking on Reference image...\n\n');
            
            % Run APPLYTOPUP
            fprintf('Correcting distortion by APPLYTOPUP...\n');
            PREPROC.dc_func_reference_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_dc_reference.nii']);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'applytopup' ...
                ' --imain=' PREPROC.func_reference_files{i} ...
                ' --inindex=' num2str(applytopup_idx) ...
                ' --topup=' PREPROC.topup.output_path ...
                ' --datain=' PREPROC.topup.dc_param_file ...
                ' --method=jac' ...
                ' --interp=spline' ...
                ' --out=' PREPROC.dc_func_reference_files{i}]);
            
            % Removing spline interpolation neg values by absolute
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'fslmaths' ...
                ' ' PREPROC.dc_func_reference_files{i} ...
                ' -abs' ...
                ' ' PREPROC.dc_func_reference_files{i}]);
            
            % Extract brain mask of reference image
            fprintf('Extracting brain of reference image and saving implicit mask...\n');
            PREPROC.dc_func_reference_files_masked{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_dc_reference_masked.nii']);
            PREPROC.dc_func_reference_files_binarymask{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_dc_reference_masked_mask.nii']);
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'bet' ...
                ' ' PREPROC.dc_func_reference_files{i} ...
                ' ' PREPROC.dc_func_reference_files_masked{i} ...
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
        ref = PREPROC.dc_func_reference_files{ref_run};
        ref_masked = PREPROC.dc_func_reference_files_masked{ref_run};
        ref_binmask = PREPROC.dc_func_reference_files_binarymask{ref_run}; 
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_dc_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        PREPROC.mean_dc_func_bold_files_masked{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '_masked.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.dc_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_dc_func_bold_files{i}]);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.mean_dc_func_bold_files{i} ...
            ' -mul ' ref_binmask ...
            ' ' PREPROC.mean_dc_func_bold_files_masked{i}]);
    
    end
    
end

% Take snapshot of distortion-corrected reference images
fprintf('Taking snapshot of distortion-corrected reference images.\n');
canlab_preproc_show_montage(PREPROC.dc_func_reference_files_masked, fullfile(PREPROC.qcdir, 'dc_func_reference_masked.png'));
drawnow;

% Take snapshot of distortion-corrected mean functional images
fprintf('Taking snapshot of distortion-corrected mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_dc_func_bold_files_masked, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold_masked.png']));
drawnow;

close all;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end