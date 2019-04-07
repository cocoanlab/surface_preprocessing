function PREPROC = humanfmri_s4_anatomical_segmentation_normalization(subject_code, study_imaging_dir, oasis_dir, varargin)

% This function does anatomical segmentation and normalization using
% antsBrainExtraction.sh and antsRegistration (ANTs).
% Note: Skip this step if CIFTIFY was done before.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s4_anatomical_segmentation_normalization(subject_code, study_imaging_dir, oasis_dir, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%   - oasis_dir          the directory containing OASIS brain template for
%                        brain extraction.
%
%
% :Optional Input:
% ::
%   - n_thread            the number of threads for ANTs (Default: 1)
%
%
% :Output:
% ::
%     PREPROC.anat_reference_file
%     PREPROC.anat_reference_file_masked
%     PREPROC.anat_reference_file_binarymask
%     PREPROC.anat_reference_file_warped
%     PREPROC.anat_reference_file_warpmat
%     PREPROC.anat_reference_file_warpfield
%     PREPROC.anat_reference_file_invwarpfield
%     PREPROC.anat_reference_file_wmseg_coreg
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
n_thread = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'n_thread'}
                n_thread = varargin{i+1};
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Anatomical segmentation and normalization', PREPROC.subject_code);
PREPROC.current_step = 's4';
PREPROC.current_step_letter = PREPROC.current_step_letter;

% Skull stripping
fprintf('Bias field correction and skull stripping...\n');
PREPROC.anat_reference_file = fullfile(PREPROC.preproc_anat_dir, 'anat_reference.nii');
PREPROC.anat_reference_file_masked = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_masked.nii');
PREPROC.anat_reference_file_binarymask = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_masked_mask.nii');
system(['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=' num2str(n_thread) ';' ...
    'antsBrainExtraction.sh' ...
    ' -k 1' ... % keep temporary files (including N4 bias corrected image)
    ' -u 0' ... % no random seeding
    ' -s nii' ...
    ' -d 3' ...
    ' -q 1' ... % float precision
    ' -a ' PREPROC.anat_files{~cellfun(@isempty, regexp(PREPROC.anat_files, 'T1w'))} ...
    ' -e ' fullfile(oasis_dir, 'T_template0.nii.gz') ...
    ' -m ' fullfile(oasis_dir, 'T_template0_BrainCerebellumProbabilityMask.nii.gz') ...
    ' -f ' fullfile(oasis_dir, 'T_template0_BrainCerebellumRegistrationMask.nii.gz') ...
    ' -o ' fullfile(PREPROC.preproc_anat_dir, 'antsBE_output_')]);
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsBE_output_N4Corrected0.nii') ' ' PREPROC.anat_reference_file]);
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsBE_output_BrainExtractionBrain.nii') ' ' PREPROC.anat_reference_file_masked]);
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsBE_output_BrainExtractionMask.nii') ' ' PREPROC.anat_reference_file_binarymask]);
system(['rm ' fullfile(PREPROC.preproc_anat_dir, 'antsBE_output_*')]);

% Normalization
PREPROC.anat_reference_file_warped = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warped.nii');
PREPROC.anat_reference_file_warpmat = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warpmat.mat');
PREPROC.anat_reference_file_warpfield = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_warpfield.nii');
PREPROC.anat_reference_file_invwarpfield = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_invwarpfield.nii');
system(['antsRegistrationSyn.sh' ...
    ' -d 3' ...
    ' -f $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz' ...
    ' -m ' PREPROC.anat_reference_file_masked ...
    ' -o ' fullfile(PREPROC.preproc_anat_dir, 'antsReg_output_') ...
    ' -n ' num2str(n_thread) ...
    ' -t s' ... % Rigid + Affine + SyN
    ' -p f' ... % float precision
    ' -j 1' ... % histogram matching
    ' -z 1']); % collapse transformation matrices
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsReg_output_0GenericAffine.mat') ' ' PREPROC.anat_reference_file_warpmat]);
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsReg_output_1Warp.nii.gz') ' ' PREPROC.anat_reference_file_warpfield '.gz']);
system(['mv ' fullfile(PREPROC.preproc_anat_dir, 'antsReg_output_1InverseWarp.nii.gz') ' ' PREPROC.anat_reference_file_invwarpfield '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_warpfield '.gz']);
system(['export FSLOUTPUTTYPE=NIFTI; fslchfiletype NIFTI ' PREPROC.anat_reference_file_invwarpfield '.gz']);
system(['rm ' fullfile(PREPROC.preproc_anat_dir, 'antsReg_output_*')]);
system(['export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=' num2str(n_thread) ';' ...
    'antsApplyTransforms' ...
    ' --verbose 1' ...
    ' --dimensionality 3' ...
    ' --float 1' ...
    ' --input ' PREPROC.anat_reference_file ...
    ' --interpolation Linear' ...
    ' --output ' PREPROC.anat_reference_file_warped ...
    ' --reference-image $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz' ...
    ' --transform ' PREPROC.anat_reference_file_warpfield ...
    ' --transform ' PREPROC.anat_reference_file_warpmat]);

% WM/CSF segmentation
fprintf('WM/CSF segmentation...\n');
system(['export FSLOUTPUTTYPE=NIFTI;' ...
    ...
    'fast' ...
    ' -o ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_fastseg') ...
    ' ' PREPROC.anat_reference_file_masked]);
PREPROC.anat_reference_file_wmseg_coreg = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_coreg.nii');
PREPROC.anat_reference_file_wmseg_nuisance = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_wmseg_nuisance.nii');
PREPROC.anat_reference_file_csfseg_nuisance = fullfile(PREPROC.preproc_anat_dir, 'anat_reference_csfseg_nuisance.nii');
system(['export FSLOUTPUTTYPE=NIFTI;' ...
    ...
    'fslmaths' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_fastseg_pve_2.nii') ...
    ' -thr 0.8' ...
    ' -bin' ...
    ' ' PREPROC.anat_reference_file_wmseg_coreg]);
system(['export FSLOUTPUTTYPE=NIFTI;' ...
    ...
    'fslmaths' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_fastseg_pve_2.nii') ...
    ' -thr 0.95' ...
    ' -bin' ...
    ' ' PREPROC.anat_reference_file_wmseg_nuisance]);
system(['export FSLOUTPUTTYPE=NIFTI;' ...
    ...
    'fslmaths' ...
    ' ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_fastseg_pve_0.nii') ...
    ' -thr 0.95' ...
    ' -bin' ...
    ' ' PREPROC.anat_reference_file_csfseg_nuisance]);
system(['rm ' fullfile(PREPROC.preproc_anat_dir, 'anat_reference_fastseg*')]);

% Erosion of WM/CSF masks for nuisance regression
fprintf('Erosion of WM/CSF masks for nuisance regression...\n');
for i = 1:5
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
    if nvox(1) < 100
        system(['rm ' PREPROC.anat_reference_file_wmseg_nuisance_erosion{i}]);
        PREPROC.anat_reference_file_wmseg_nuisance_erosion(i) = [];
        break;
    end
end
for i = 1:2
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
    if nvox(1) < 100
        system(['rm ' PREPROC.anat_reference_file_csfseg_nuisance_erosion{i}]);
        PREPROC.anat_reference_file_csfseg_nuisance_erosion(i) = [];
        break;
    end
end

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end
