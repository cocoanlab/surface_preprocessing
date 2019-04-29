function PREPROC = humanfmri_s10_denoising(subject_code, study_imaging_dir, varargin)

% This function de-noises linear trend, WM/CSF, and extremely high or low
% frequency signals using 3dTproject function (AFNI).
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s10_denoising(subject_code, study_imaging_dir, varargin)
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
%
%
% :Output:
% ::
%     PREPROC.n_func_bold_files
%     PREPROC.mean_n_func_bold_files
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
detrend_dim = 0;
bandpass = [0 9999];
tr = {[]};
do_movement = false;
do_wmcsf = false;
do_custom = false;
% Add 'censor'
do_addmean = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'detrend'}
                detrend_dim = varargin{i+1};
            case {'bandpass'}
                bandpass = varargin{i+1};
            case {'tr'}
                tr = varargin{i+1};
                if ~iscell(tr)
                    tr = {tr};
                end
            case {'movement'}
                do_movement = true;
                movement_method = varargin{i+1};
            case {'wmcsf'}
                do_wmcsf = true;
                wmcsf_method = varargin{i+1};
            case {'custom_nuisance'}
                do_custom = true;
                custom_nuisance = varargin{i+1};
            case {'addmean'}
                do_addmean = true;
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Denoising', PREPROC.subject_code);
PREPROC.current_step = 's10';
PREPROC.current_step_letter = ['n' PREPROC.current_step_letter];

if numel(tr) == 1
    tr = repmat(tr, numel(PREPROC.func_bold_files), 1);
end
tr = tr(:);

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
        
        % Setting
        PREPROC.n_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        
        nuisance_mat = [];
        
        if detrend_dim ~= 0
            fprintf('*** Polynomial trend: Up to degree %d ***\n', detrend_dim);
        end
        
        if do_movement
            mvmt_dat = importdata(PREPROC.mvmt_param_files_demeaned{i});
            switch movement_method
                case '6p'
                    nuisance_mat = [nuisance_mat mvmt_dat];
                case '24p'
                    mvmt_dat_diff = [zeros(1,size(mvmt_dat,2)); diff(mvmt_dat)];
                    nuisance_mat = [nuisance_mat mvmt_dat mvmt_dat_diff mvmt_dat.^2 mvmt_dat_diff.^2];
            end
        end
        
        if do_wmcsf
            if regexp(PREPROC.current_step_letter, 'w')
                wm_mask = which('canonical_white_matter.img');
                csf_mask = which('canonical_ventricles.img');
            else
                wm_mask = PREPROC.coregistered_wmseg_nuisance_ero;
                csf_mask = PREPROC.coregistered_csfseg_nuisance_ero;
            end
            [~, temp_denoising_dir] = system('mktemp -d');
            temp_denoising_dir = strtrim(temp_denoising_dir);
            WM_roi = fullfile(temp_denoising_dir, 'WM_roi.1D');
            CSF_roi = fullfile(temp_denoising_dir, 'CSF_roi.1D');
            system(['3dmaskdump' ...
                ' -noijk' ...
                ' -mask ' wm_mask ...
                ' ' PREPROC.i_func_bold_files{i} ...
                ' > ' WM_roi]);
            system(['3dmaskdump' ...
                ' -noijk' ...
                ' -mask ' csf_mask ...
                ' ' PREPROC.i_func_bold_files{i} ...
                ' > ' CSF_roi]);
            WM_dat = importdata(WM_roi);
            CSF_dat = importdata(CSF_roi);
            system(['rm -r ' temp_denoising_dir]);
                
            switch wmcsf_method
                case 'mean'
                    fprintf('*** WM/CSF: Mean signal ***\n');
                    nuisance_mat = [nuisance_mat mean(WM_dat)' mean(CSF_dat)'];
                case 'acompcor'
                    fprintf('*** WM/CSF: aCompcor ***\n');
                    [~, WM_PC] = pca(WM_dat');
                    [~, CSF_PC] = pca(CSF_dat');
                    nuisance_mat = [nuisance_mat WM_PC(:,1:5) CSF_PC(:,1:5)];
                case 'acompcor50'
                    fprintf('*** WM/CSF: aCompcor50 ***\n');
                    [~, WM_PC, ~, ~, WM_PC_explained] = pca(WM_dat');
                    WM_PC50_idx = find(cumsum(WM_PC_explained) >= 50, 1);
                    [~, CSF_PC, ~, ~, CSF_PC_explained] = pca(CSF_dat');
                    CSF_PC50_idx = find(cumsum(CSF_PC_explained) >= 50, 1);
                    fprintf('WM: PC%.2d ~ PC%.2d explains %.4f percent of the variance.\n', ...
                        1, WM_PC50_idx, sum(WM_PC_explained(1:WM_PC50_idx)));
                    fprintf('CSF: PC%.2d ~ PC%.2d explains %.4f percent of the variance.\n', ...
                        1, CSF_PC50_idx, sum(CSF_PC_explained(1:CSF_PC50_idx)));
                    nuisance_mat = [nuisance_mat WM_PC(:,1:WM_PC50_idx) CSF_PC(:,1:CSF_PC50_idx)];
            end
            WM_dat = [];
            CSF_dat = [];
        end
        
        if do_custom
            fprintf('*** Custom nuisance variables ***\n');
            nuisance_mat = [nuisance_mat custom_nuisance{i}];
        end
        
        PREPROC.nuisance_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['nuisance_' b '.1D']);
        dlmwrite(PREPROC.nuisance_files{i}, nuisance_mat, 'delimiter', '\t');
        
        if isempty(tr{i})
            load(PREPROC.func_bold_dicomheader_files{i});
            tr{i} = h.RepetitionTime; % msec
        end
        
        system(['3dTproject' ...
            ' -verb' ...
            ' -input ' PREPROC.i_func_bold_files{i} ...
            ' -mask ' func_ref_mask ...
            ' -ort ' PREPROC.nuisance_files{i} ...
            ' -polort ' num2str(detrend_dim) ...
            ' -bandpass ' num2str(bandpass(1)) ' ' num2str(bandpass(2)) ...
            ' -TR ' num2str(tr{i} / 1000) ...
            ' -prefix ' PREPROC.n_func_bold_files{i}]);
        
        if do_addmean
            fprintf('*** Adding mean back to the denoised data ***\n');
            system(['export FSLOUTPUTTYPE=NIFTI;' ...
                ...
                'fslmaths' ...
                ' ' PREPROC.n_func_bold_files{i} ...
                ' -add' ...
                ' ' PREPROC.mean_i_func_bold_files{i} ...
                ' ' PREPROC.n_func_bold_files{i}]);
        end
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_n_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.n_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_n_func_bold_files{i}]);
        
    end
    
end

% Take snapshot of denoised mean functional images
fprintf('Taking snapshot of denoised mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_n_func_bold_files, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold.png']));
drawnow;

close all;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end