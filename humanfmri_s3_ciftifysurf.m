function PREPROC = humanfmri_s3_ciftifysurf(subject_code, study_imaging_dir, ciftify_basedir, varargin)

% This function converts surface data that is reconstructed by Freesurfer
% ('recon-all' function') to CIFTI format, using ciftify_recon_all
% function in CIFTIFY toolbox ('https://github.com/edickie/ciftify').
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s3_ciftifysurf(subject_code, study_imaging_dir, ciftify_basedir)
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
%   - n_ero_limit        limit of number of erosion.
%                        (e.g., [3 1]: 3 and 1 times for WM and CSF)
%                        (default: [5 2]).
%   - nvox_ero_limit     limit of number of voxels for erosion
%                        (default: 100).
%
%
% :Output:
% ::
%     PREPROC.study_cifti_dir
%     PREPROC.anat_reference_file
%     PREPROC.anat_reference_file_masked
%     PREPROC.anat_reference_file_binarymask
%     PREPROC.anat_reference_file_warped
%     PREPROC.anat_reference_file_warped_binarymask
%     PREPROC.anat_reference_file_warpfield
%     PREPROC.anat_reference_file_invwarpfield
%     PREPROC.anat_reference_file_wmseg_for_coreg
%     PREPROC.anat_reference_file_wmseg_nuisance
%     PREPROC.anat_reference_file_csfseg_nuisance
%     PREPROC.anat_reference_file_wmseg_nuisance_erosion
%     PREPROC.anat_reference_file_csfseg_nuisance_erosion
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
nvox_ero_limit = 100;
n_ero_limit = [5 2];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'n_ero_limit'}
                n_ero_limit = varargin{i+1};
            case {'nvox_ero_limit'}
                nvox_ero_limit = varargin{i+1};
        end
    end
end

PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Converting surface to CIFTI format', PREPROC.subject_code);
PREPROC.current_step = 's3';
PREPROC.current_step_letter = PREPROC.current_step_letter;

PREPROC.anat_normalization_method = 'FSL';

PREPROC.cifti_dir = fullfile(PREPROC.study_imaging_dir, 'cifti');
if ~exist(PREPROC.cifti_dir, 'dir')
    mkdir(PREPROC.cifti_dir);
end

% CIFTIFY work
system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
    'export PATH=' ciftify_basedir '/msm:' ciftify_basedir '/ciftify/bin:$PATH;' ...
    'export PYTHONPATH=' ciftify_basedir ':$PYTHONPATH;' ...
    'export CIFTIFY_TEMPLATES=' ciftify_basedir '/data;' ...
    'export CIFTIFY_WORKDIR=' PREPROC.cifti_dir ';' ...
    ...
    'ciftify_recon_all' ...
    ' --verbose' ...
    ' ' PREPROC.subject_code]);

% CIFTIFY visualization
system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
    'export PATH=' ciftify_basedir '/msm:' ciftify_basedir '/ciftify/bin:$PATH;' ...
    'export PYTHONPATH=' ciftify_basedir ':$PYTHONPATH;' ...
    'export CIFTIFY_TEMPLATES=' ciftify_basedir '/data;' ...
    'export CIFTIFY_WORKDIR=' PREPROC.cifti_dir ';' ...
    ...
    'cifti_vis_recon_all' ...
    ' subject' ...
    ' --verbose' ...
    ' ' PREPROC.subject_code]);



% Copy anatomical reference files from CIFTIFY
fprintf('Copying anatomical reference files from CIFTIFY.\n');
PREPROC.anat_reference_file = fullfile(PREPROC.preproc_anat_dir, 'anat_reference.nii');
PREPROC.anat_reference_file_masked = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_masked.nii');
PREPROC.anat_reference_file_binarymask = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_masked_mask.nii');
PREPROC.anat_reference_file_warped = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warped.nii');
PREPROC.anat_reference_file_warped_binarymask = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warped_masked_mask.nii');
PREPROC.anat_reference_file_warpfield = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warpfield.nii');
PREPROC.anat_reference_file_invwarpfield = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_invwarpfield.nii');
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'T1w.nii.gz') ' ' PREPROC.anat_reference_file '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'T1w_brain.nii.gz') ' ' PREPROC.anat_reference_file_masked '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'brainmask_fs.nii.gz') ' ' PREPROC.anat_reference_file_binarymask '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'MNINonLinear', 'T1w.nii.gz') ' ' PREPROC.anat_reference_file_warped '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'MNINonLinear', 'brainmask_fs.nii.gz') ' ' PREPROC.anat_reference_file_warped_binarymask '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'MNINonLinear', 'xfms', 'T1w2Standard_warp_noaffine.nii.gz') ' ' PREPROC.anat_reference_file_warpfield '.gz']);
system(['cp ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'MNINonLinear', 'xfms', 'Standard2T1w_warp_noaffine.nii.gz') ' ' PREPROC.anat_reference_file_invwarpfield '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_masked '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_binarymask '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_warped '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_warped_binarymask '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_warpfield '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_invwarpfield '.gz']);

% Generate GM segmentation using wmparc output from CIFTIFY
PREPROC.anat_reference_file_gmseg = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_gmseg.nii');
system(['mri_binarize' ...
    ' --i ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'aparc+aseg.nii.gz') ...
    ' --gm' ...
    ' --o ' PREPROC.anat_reference_file_gmseg]);

% Generate WM segmentation using wmparc output from CIFTIFY
% ''' use the wmparc file in the anat folder to define the wm mask
%     will do so by combining
%     LEFT-CEREBRAL-WHITE-MATTER
%     2 245 245 245 255
%     LEFT-CEREBELLUM-WHITE-MATTER
%     7 220 248 164 255
%     RIGHT-CEREBRAL-WHITE-MATTER
%     41 0 225 0 255
%     RIGHT-CEREBELLUM-WHITE-MATTER
%     46 220 248 164 255
%     CC_POSTERIOR
%     251 0 0 64 255
%     CC_MID_POSTERIOR
%     252 0 0 112 255
%     CC_CENTRAL
%     253 0 0 160 255
%     CC_MID_ANTERIOR
%     254 0 0 208 255
%     CC_ANTERIOR
%     255 0 0 255 255
%     all the 3000*s and 4000*s 5001 5002
%     but there is also the question of the deep gray matter (that can look white?)
%     BRAINSTEM, PALLIDUM, THALAMUS (thalamus is half and half), VENTRALDC
%     LEFT-PALLIDUM
%     13 12 48 255 255
%     LEFT-VENTRALDC
%     28 165 42 42 255
%     BRAIN-STEM
%     16 119 159 176 255
%     RIGHT-PALLIDUM
%     52 13 48 255 255
%     RIGHT-VENTRALDC
%     60 165 42 42 255
%     '''
PREPROC.anat_reference_file_wmseg_for_coreg = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_for_coreg.nii');
PREPROC.anat_reference_file_wmseg_nuisance = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_nuisance.nii');
system(['wb_command' ...
    ' -volume-math "(x == 2 || x == 7 || x == 41 || x == 46)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part1.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "(x == 251 || x == 252 || x == 253 || x == 254 || x == 255)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part2.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "(x > 2999 && x < 5005)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part3.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "(x == 13 || x == 28 || x == 16 || x == 52 || x == 60)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part4.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "((a + b + c + d) > 0)"' ...
    ' ' PREPROC.anat_reference_file_wmseg_for_coreg ...
    ' -var a ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part1.nii') ...
    ' -var b ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part2.nii') ...
    ' -var c ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part3.nii') ...
    ' -var d ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part4.nii')]);
system(['wb_command' ...
    ' -volume-math "((a + b + c) > 0)"' ...
    ' ' PREPROC.anat_reference_file_wmseg_nuisance ...
    ' -var a ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part1.nii') ...
    ' -var b ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part2.nii') ...
    ' -var c ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part3.nii')]);
system(['rm ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_part*.nii')]);
for i = 1:n_ero_limit(1)
    PREPROC.anat_reference_file_wmseg_nuisance_erosion{i, 1} = fullfile(PREPROC.preproc_anat_dir, ['anat_reference_wmseg_nuisance_ero' num2str(i) '.nii']);
    if i == 1
        input_dat = PREPROC.anat_reference_file_wmseg_nuisance;
    else
        input_dat = PREPROC.anat_reference_file_wmseg_nuisance_erosion{i-1};
    end
    output_dat = PREPROC.anat_reference_file_wmseg_nuisance_erosion{i};
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'fslmaths ' input_dat ' -ero ' output_dat]);
    [~, nvox] = system(['fslstats ' output_dat ' -V']);
    nvox = str2num(nvox);
    if nvox(1) < nvox_ero_limit
        system(['rm ' PREPROC.anat_reference_file_wmseg_nuisance_erosion{i}]);
        PREPROC.anat_reference_file_wmseg_nuisance_erosion(i) = [];
        break;
    end
end


% Generate CSF segmentation using wmparc output from CIFTIFY
% ''' use the wmparc file in the anat folder to define the csf mask
%     will do so by combining
%     LEFT-LATERAL-VENTRICLE
%     4 120 18 134 0
%     LEFT-INF-LAT-VENT
%     5 196 58 250 0
%     3RD-VENTRICLE
%     14 204 182 142 0
%     4TH-VENTRICLE         
%     15 42 204 164 0
%     CSF                                     
%     24 60 60 60 0
%     LEFT-CHOROID-PLEXUS
%     31 0 200 200 0
%     RIGHT-LATERAL-VENTRICLE           
%     43 120 18 134 0
%     RIGHT-INF-LAT-VENT                  
%     44 196 58 250 0
%     RIGHT-CHOROID-PLEXUS                
%     63 0 200 221 0
%     5TH-VENTRICLE
%     72 120 190 150 0
%     '''
PREPROC.anat_reference_file_csfseg_nuisance = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_nuisance.nii');
system(['wb_command' ...
    ' -volume-math "(x == 4 || x == 5 || x == 43 || x == 44 || x == 31 || x == 63)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_part1.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "(x == 14 || x == 15 || x == 72 || x == 24)"' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_part2.nii') ...
    ' -var x ' fullfile(PREPROC.cifti_dir, PREPROC.subject_code, 'T1w', 'wmparc.nii.gz')]);
system(['wb_command' ...
    ' -volume-math "((a + b) > 0)"' ...
    ' ' PREPROC.anat_reference_file_csfseg_nuisance ...
    ' -var a ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_part1.nii') ...
    ' -var b ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_part2.nii')]);
system(['rm ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_part*.nii')]);
for i = 1:n_ero_limit(2)
    PREPROC.anat_reference_file_csfseg_nuisance_erosion{i, 1} = fullfile(PREPROC.preproc_anat_dir, ['anat_reference_csfseg_nuisance_ero' num2str(i) '.nii']);
    if i == 1
        input_dat = PREPROC.anat_reference_file_csfseg_nuisance;
    else
        input_dat = PREPROC.anat_reference_file_csfseg_nuisance_erosion{i-1};
    end
    output_dat = PREPROC.anat_reference_file_csfseg_nuisance_erosion{i};
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'fslmaths ' input_dat ' -ero ' output_dat]);
    [~, nvox] = system(['fslstats ' output_dat ' -V']);
    nvox = str2num(nvox);
    if nvox(1) < nvox_ero_limit
        system(['rm ' PREPROC.anat_reference_file_csfseg_nuisance_erosion{i}]);
        PREPROC.anat_reference_file_csfseg_nuisance_erosion(i) = [];
        break;
    end
end

save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end
