function PREPROC = humanfmri_s13_ciftifyfmri(subject_code, study_imaging_dir, ciftify_basedir, varargin)

% This function converts functional data to CIFTI format, using ciftify_recon_all
% function in CIFTIFY toolbox ('https://github.com/edickie/ciftify').
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s13_ciftifyfmri(subject_code, study_imaging_dir, ciftify_basedir, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%   - ciftify_basedir    the cifitify toolbox directory.
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
%     PREPROC.cifti_func_dir
%     PREPROC.cifti_func_bold_files
%     PREPROC.cifti_anat_midthickness_L_file
%     PREPROC.cifti_anat_midthickness_R_file
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
print_header('Converting functional data to CIFTI format', PREPROC.subject_code);
PREPROC.current_step = 's13';
PREPROC.current_step_letter = PREPROC.current_step_letter;

if ~strcmp(PREPROC.anat_normalization_method, 'FSL')
    error('Please run step ''s3'' before you run this function.')
end

for i = 1:numel(PREPROC.func_bold_files)
    
    if ~do_select_run || ismember(i, run_num)
        
        fprintf('\n\nWorking on Run %d...\n\n', i);
        [~, b] = fileparts(PREPROC.func_bold_files{i});
        
        if PREPROC.ref_first_run
            ref_run = 1;
        else
            ref_run = i;
        end
        ref = PREPROC.w_func_reference_file{ref_run};
        
        % CIFTIFY work
        [~, functype] = fileparts(PREPROC.dicom_func_bold_dir{i});
        system(['export FSLOUTPUTTYPE=NIFTI_GZ;' ...
            ...
            'export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            'export PATH=' ciftify_basedir '/msm:' ciftify_basedir '/ciftify/bin:$PATH;' ...
            'export PYTHONPATH=' ciftify_basedir ':$PYTHONPATH;' ...
            'export CIFTIFY_TEMPLATES=' ciftify_basedir '/data;' ...
            'export CIFTIFY_WORKDIR=' PREPROC.cifti_dir ';' ...
            ...
            'ciftify_subject_fmri' ...
            ' --verbose' ...
            ' --func-ref ' ref ...
            ' --already-in-MNI' ...
            ' --OutputSurfDiagnostics' ...
            ' --SmoothingFWHM ' num2str(fwhm) ... % 5mm FWHM
            ' ' PREPROC.w_func_bold_files{i} ...
            ' ' PREPROC.subject_code ...
            ' ' functype]);
        
        % CIFTIFY visualization
        system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            'export PATH=' ciftify_basedir '/msm:' ciftify_basedir '/ciftify/bin:$PATH;' ...
            'export PYTHONPATH=' ciftify_basedir ':$PYTHONPATH;' ...
            'export CIFTIFY_TEMPLATES=' ciftify_basedir '/data;' ...
            'export CIFTIFY_WORKDIR=' PREPROC.cifti_dir ';' ...
            ...
            'cifti_vis_fmri' ...
            ' subject' ...
            ' --verbose' ...
            ' --SmoothingFWHM ' num2str(fwhm) ... % 5mm FWHM
            ' --meanfunc ' PREPROC.mean_w_func_bold_files{i} ...
            ' ' functype ...
            ' ' PREPROC.subject_code]);
        
        PREPROC.cifti_func_dir{i, 1} = fullfile(PREPROC.cifti_dir, subject_code, 'MNINonLinear', 'Results', functype);
        PREPROC.cifti_func_bold_files{i, 1} = fullfile(PREPROC.cifti_func_dir{i}, [functype '_Atlas_s5.dtseries.nii']);
        
    end
    
end

PREPROC.cifti_anat_midthickness_L_file = fullfile(PREPROC.cifti_dir, subject_code, 'MNINonLinear', 'fsaverage_LR32k', [subject_code '.L.midthickness.32k_fs_LR.surf.gii']);
PREPROC.cifti_anat_midthickness_R_file = fullfile(PREPROC.cifti_dir, subject_code, 'MNINonLinear', 'fsaverage_LR32k', [subject_code '.R.midthickness.32k_fs_LR.surf.gii']);

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end
