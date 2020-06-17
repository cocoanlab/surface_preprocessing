function PREPROC = humanfmri_s6_motion_correction(subject_code, study_imaging_dir, varargin)

% This function does motion correction (realignment) on functional data
% using 3dvolreg function (AFNI). If slice timing correction is performed,
% movement parameter is not calculated (since it should be calculated
% before slice timing correction).
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s6_motion_correction(subject_code, study_imaging_dir, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001')
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging')
%
%
% :Optional Input:
% ::
%   - run_num            runs to include. ex) [1 2 4 5], ...
%
% :Output:
% ::
%     PREPROC.r_func_bold_files
%     PREPROC.mvmt_param_files
%     PREPROC.mvmt_matrix_files
%     PREPROC.nuisance
%     PREPROC.mean_r_func_bold_files
%     PREPROC.mean_r_func_bold_files_masked
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

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Motion correction (realignment)', PREPROC.subject_code);
PREPROC.current_step = 's6';
PREPROC.current_step_letter = ['r' PREPROC.current_step_letter];

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        if PREPROC.ref_first_run
            ref_run = 1;
        else
            ref_run = i;
        end
        ref = PREPROC.func_reference_files{ref_run};
        ref_masked = PREPROC.func_reference_files_masked{ref_run};
        ref_binmask = PREPROC.func_reference_files_binarymask{ref_run};
        
        % Motion correction (realignment)
        PREPROC.r_func_bold_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [PREPROC.current_step_letter b '.nii']);
        if strcmp(PREPROC.current_step_letter, 'ra')
            fprintf('Slice timing-corrected images.\n');
            fprintf('Motion correction without calculating movement parameters...\n');
            system(['3dvolreg' ...
                ' -base ' ref ...
                ' -verbose' ...
                ' -Fourier' ...
                ' -twopass' ...
                ' -zpad 4' ...
                ' -float' ...
                ' -prefix ' PREPROC.r_func_bold_files{i} ...
                ' ' PREPROC.a_func_bold_files{i}]);
        elseif strcmp(PREPROC.current_step_letter, 'r')
            fprintf('Motion correction with calculating movement parameters...\n');
            PREPROC.mvmt_param_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['rp_' b '.1D']);
            PREPROC.mvmt_matrix_files{i, 1} = fullfile(PREPROC.preproc_func_dir, ['rp_mat_' b '.1D']);
            system(['3dvolreg' ...
                ' -base ' ref ...
                ' -verbose' ...
                ' -Fourier' ...
                ' -twopass' ...
                ' -zpad 4' ...
                ' -1Dfile ' PREPROC.mvmt_param_files{i} ... % Roll/Pitch/Yaw (deg, counter-clockwise) & dS/dL/dP (mm)
                ' -1Dmatrix_save ' PREPROC.mvmt_matrix_files{i} ...
                ' -float' ...
                ' -prefix ' PREPROC.r_func_bold_files{i} ...
                ' ' PREPROC.func_bold_files{i}]);
        end
        
        mvmtdata = importdata(PREPROC.mvmt_param_files{i});
        if isstruct(mvmtdata)
            mvmtdata = mvmtdata.data;
        end
        PREPROC.mvmt_param_files_demeaned{i, 1} = fullfile(PREPROC.preproc_func_dir, ['rp_demeaned_' b '.1D']);
        dlmwrite(PREPROC.mvmt_param_files_demeaned{i}, scale(mvmtdata, 1), 'delimiter', '\t');
        
        % save plot
        fprintf('Taking snapshot of movement parameters.\n');
        create_figure('mvmt', 2, 1);
        
        subplot(2,1,1);
        plot([-mvmtdata(:,5), -mvmtdata(:,6), mvmtdata(:,4)]);
        legend('x', 'y', 'z');
        xlabel('volume');
        ylabel(sprintf('displacement\n(mm)'));
        
        subplot(2,1,2);
        plot([mvmtdata(:,2), mvmtdata(:,1), mvmtdata(:,3)]);
        legend('pitch', 'roll', 'yaw');
        xlabel('volume');
        ylabel(sprintf('displacement\n(deg, counter-clockwise)'));
        
        sz = get(0, 'screensize'); % Wani added two lines to make this visible (but it depends on the size of the monitor)
        set(gcf, 'Position', [sz(3)*.02 sz(4)*.05 sz(3) *.45 sz(4)*.85]);
        drawnow;
        
        mvmt_qcfile = fullfile(PREPROC.qcdir, ['qc_mvmt_' b '.png']); % Scott added some lines to actually save the spike images
        saveas(gcf, mvmt_qcfile);
        close all;
        
        % Save mean functional image
        fprintf('Saving mean functional image...\n');
        PREPROC.mean_r_func_bold_files{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '.nii']);
        PREPROC.mean_r_func_bold_files_masked{i, 1} = fullfile(PREPROC.preproc_mean_func_dir, ['mean_' PREPROC.current_step_letter b '_masked.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.r_func_bold_files{i} ...
            ' -Tmean' ...
            ' ' PREPROC.mean_r_func_bold_files{i}]);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'fslmaths' ...
            ' ' PREPROC.mean_r_func_bold_files{i} ...
            ' -mul ' ref_binmask ...
            ' ' PREPROC.mean_r_func_bold_files_masked{i}]);
        
    end
    
end

% Take snapshot of motion-corrected mean functional images
fprintf('Taking snapshot of motion-corrected mean functional images.\n');
canlab_preproc_show_montage(PREPROC.mean_r_func_bold_files_masked, fullfile(PREPROC.qcdir, ['mean_' PREPROC.current_step_letter '_func_bold_masked.png']));
drawnow;

close all;

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end
