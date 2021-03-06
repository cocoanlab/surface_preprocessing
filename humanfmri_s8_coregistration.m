function PREPROC = humanfmri_s8_coregistration(subject_code, study_imaging_dir, varargin)

% This function does coregistration between anatomical T1 image and
% functional reference image using FLIRT BBR (FSL, for pre-alignment).
% If freesurfer surface reconstruction was done before, then bbregister
% (freesurfer) was additionally used for refinement for coregistration.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_s8_coregistration(subject_code, study_imaging_dir, varargin)
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
%   - no_check_reg       no check regisration. if you want to run all the
%                        subject without any interaction, this will be helpful.
%
%
% :Output:
% ::
%     PREPROC.transform_EPI_to_T1
%     PREPROC.transform_T1_to_EPI
%     (optional) PREPROC.transform_EPI_to_T1_ITK
%     (optional) PREPROC.transform_T1_to_EPI_ITK
%     PREPROC.coregistered_func_reference_files
%     PREPROC.coregistered_func_reference_files_masked
%     PREPROC.coregistered_func_reference_files_binarymask
%     PREPROC.coregistered_wmseg
%     PREPROC.coregistered_csfseg
%     saves qc_images/coreg_func_reference_masked.png 
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
do_check = true;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'no_check_reg'}
                do_check = false;
        end
    end
end


PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
print_header('Coregistration', PREPROC.subject_code);
PREPROC.current_step = 's8';
PREPROC.current_step_letter = PREPROC.current_step_letter;

for i = 1:numel(PREPROC.func_reference_files)
    
    fprintf('\n\nWorking on Reference %d...\n\n', i);
    [~, b] = fileparts(PREPROC.func_bold_files{i});
    
    % Step 1: Register EPI to T1 using FSL (BBR cost function)
    fprintf('Step 1: Registering EPI to T1 using FSL (BBR cost function)...\n');
    
    % Pre-alignment
    fprintf('Pre-alignment...\n');
    if regexp(PREPROC.current_step_letter, 'dc')
        ref = PREPROC.dc_func_reference_files{i};
    else
        if regexp(PREPROC.current_step_letter, 'r')
            ref = PREPROC.func_reference_files{i};
        end
    end
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'flirt' ...
        ' -in ' ref ...
        ' -ref ' PREPROC.anat_reference_file_masked ...
        ' -dof 6' ...
        ' -omat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslinit.mat']) ...
        ' -out ' fullfile(PREPROC.preproc_func_dir, [b '_coregistered_EPI_to_T1_fslinit.nii'])]);
    % FSL BBR
    fprintf('Running BBR...\n');
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'flirt' ...
        ' -in ' ref ...
        ' -ref ' PREPROC.anat_reference_file ...
        ' -dof 6' ...
        ' -cost bbr' ...
        ' -wmseg ' PREPROC.anat_reference_file_wmseg_for_coreg ...
        ' -init ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslinit.mat']) ...
        ' -omat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr.mat']) ...
        ' -schedule $FSLDIR/etc/flirtsch/bbr.sch' ...
        ' -out ' fullfile(PREPROC.preproc_func_dir, [b '_coregistered_EPI_to_T1_fslbbr.nii'])]);
    
    if strcmp(PREPROC.anat_normalization_method, 'FSL')
        
        % Step 2: Refine registration of EPI to T1 using freesurfer (bbregister)
        fprintf('Step 2: Refining registration of EPI to T1 using freesurfer (bbregister)...\n');
        
        % Calculate difference of orientation between FSL and freesurfer: FSL -> freesurfer
        system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            ...
            'tkregister2' ...
            ' --s ' PREPROC.subject_code ...
            ' --noedit' ...
            ' --mov ' PREPROC.anat_reference_file ...
            ' --targ ' fullfile(PREPROC.surfrecon_dir, PREPROC.subject_code, 'mri', 'T1.mgz') ...
            ' --reg ' fullfile(PREPROC.preproc_func_dir, [b '_transform_fsl_to_freesurfer.dat']) ...
            ' --regheader' ...
            ' --fslregout ' fullfile(PREPROC.preproc_func_dir, [b '_transform_fsl_to_freesurfer.mat'])]);
        
        %     transform_fsl_to_freesurfer =
        %
        %         1.0000         0         0         0
        %              0         0   -1.0000  255.0000
        %              0    1.0000         0         0
        %              0         0         0    1.0000
        
        % Calculate difference of orientation between FSL and freesurfer: freesurfer -> FSL
        system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            ...
            'tkregister2' ...
            ' --s ' PREPROC.subject_code ...
            ' --noedit' ...
            ' --mov ' fullfile(PREPROC.surfrecon_dir, PREPROC.subject_code, 'mri', 'T1.mgz') ...
            ' --targ ' PREPROC.anat_reference_file ...
            ' --reg ' fullfile(PREPROC.preproc_func_dir, [b '_transform_freesurfer_to_fsl.dat']) ...
            ' --regheader' ...
            ' --fslregout ' fullfile(PREPROC.preproc_func_dir, [b '_transform_freesurfer_to_fsl.mat'])]);
        
        %     transform_freesurfer_to_fsl =
        %
        %         1.0000         0         0         0
        %              0         0    1.0000   -0.0000
        %              0   -1.0000         0  255.0000
        %              0         0         0    1.0000
        
        % Generate initial coregistration matrix (FSL BBR) in freesurfer coordinates and TkRegister format
        system(['convert_xfm' ...
            ' ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr.mat']) ... % EPI to T1 with FSL BBR
            ' -concat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_fsl_to_freesurfer.mat']) ... % FSL to freesurfer
            ' -omat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr_freesurfer.mat'])]); % EPI to T1 with FSL BBR in freesurfer coordinates
        system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            ...
            'tkregister2' ...
            ' --s ' PREPROC.subject_code ...
            ' --noedit' ...
            ' --mov ' ref ...
            ' --targ ' fullfile(PREPROC.surfrecon_dir, PREPROC.subject_code, 'mri', 'T1.mgz') ...
            ' --fsl ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr_freesurfer.mat']) ...
            ' --reg ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr_freesurfer.dat'])]);
        % bbregister
        system(['export SUBJECTS_DIR=' PREPROC.surfrecon_dir ';' ...
            ...
            'bbregister' ...
            ' --s ' PREPROC.subject_code ...
            ' --mov ' ref ...
            ' --bold' ...
            ' --init-reg ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr_freesurfer.dat']) ...
            ' --reg ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_bbregister.dat']) ...
            ' --fslmat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_bbregister.mat']) ...
            ' --o ' fullfile(PREPROC.preproc_func_dir, [b '_coregistered_EPI_to_T1_bbregister.nii'])]);
        % Concatenating transformation matricies for coregistration
        PREPROC.transform_EPI_to_T1{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1.mat']);
        system(['convert_xfm' ...
            ' ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_bbregister.mat']) ... % EPI to T1 with bbregister
            ' -concat ' fullfile(PREPROC.preproc_func_dir, [b '_transform_freesurfer_to_fsl.mat']) ... % freesurfer to FSL
            ' -omat ' PREPROC.transform_EPI_to_T1{i}]); % EPI to T1 coregistration matrix!!!
        PREPROC.transform_T1_to_EPI{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_T1_to_EPI.mat']);
        system(['convert_xfm' ...
            ' ' PREPROC.transform_EPI_to_T1{i} ...
            ' -inverse' ...
            ' -omat ' PREPROC.transform_T1_to_EPI{i}]); % T1 to EPI coregistration matrix!!!
        
    elseif strcmp(PREPROC.anat_normalization_method, 'ANTS')
        
        PREPROC.transform_EPI_to_T1{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1.mat']);
        system(['cp ' fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_fslbbr.mat']) ' ' PREPROC.transform_EPI_to_T1{i}]);
        PREPROC.transform_T1_to_EPI{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_T1_to_EPI.mat']);
        system(['convert_xfm' ...
            ' ' PREPROC.transform_EPI_to_T1{i} ...
            ' -inverse' ...
            ' -omat ' PREPROC.transform_T1_to_EPI{i}]); % T1 to EPI coregistration matrix!!!
        
        PREPROC.transform_EPI_to_T1_ITK{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_EPI_to_T1_ITK.txt']);
        system(['wb_command' ...
            ' -convert-affine' ...
            ' -from-flirt' ...
            ' ' PREPROC.transform_EPI_to_T1{i} ...
            ' ' ref ...
            ' ' PREPROC.anat_reference_file ...
            ' -to-itk ' PREPROC.transform_EPI_to_T1_ITK{i}]);
        fid = fopen(PREPROC.transform_EPI_to_T1_ITK{i}, 'rt');
        ftext = char(fread(fid)');
        fclose(fid);
        fid = fopen(PREPROC.transform_EPI_to_T1_ITK{i}, 'wt');
        fprintf(fid, strrep(ftext, 'MatrixOffsetTransformBase', 'AffineTransform'));
        fclose(fid);
        PREPROC.transform_T1_to_EPI_ITK{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_transform_T1_to_EPI_ITK.txt']);
        system(['wb_command' ...
            ' -convert-affine' ...
            ' -from-flirt' ...
            ' ' PREPROC.transform_T1_to_EPI{i} ...
            ' ' PREPROC.anat_reference_file ...
            ' ' ref ...
            ' -to-itk ' PREPROC.transform_T1_to_EPI_ITK{i}]);
        fid = fopen(PREPROC.transform_T1_to_EPI_ITK{i}, 'rt');
        ftext = char(fread(fid)');
        fclose(fid);
        fid = fopen(PREPROC.transform_T1_to_EPI_ITK{i}, 'wt');
        fprintf(fid, strrep(ftext, 'MatrixOffsetTransformBase', 'AffineTransform'));
        fclose(fid);
        
    end
    
    % Generate coreigstered functional reference image
    PREPROC.coregistered_func_reference_files{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_coreg_reference.nii']);
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'flirt' ...
        ' -in ' ref ...
        ' -ref ' PREPROC.anat_reference_file ...
        ' -applyxfm -init ' PREPROC.transform_EPI_to_T1{i} ...
        ' -out ' PREPROC.coregistered_func_reference_files{i}]);
    % Extract brain mask of reference image
    fprintf('Extracting brain of reference image and saving implicit mask...\n');
    PREPROC.coregistered_func_reference_files_masked{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_coreg_reference_masked.nii']);
    PREPROC.coregistered_func_reference_files_binarymask{i, 1} = fullfile(PREPROC.preproc_func_dir, [b '_coreg_reference_masked_mask.nii']);
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'bet' ...
        ' ' PREPROC.coregistered_func_reference_files{i} ...
        ' ' PREPROC.coregistered_func_reference_files_masked{i} ...
        ' -f 0.3' ... % ICA-AROMA recommendation
        ' -n' ...
        ' -m' ...
        ' -R']);
    % Generate coreigstered anatomical segmentation data
    if isfield(PREPROC, 'anat_reference_file_gmseg')
        PREPROC.coregistered_gmseg{i, 1} = fullfile(PREPROC.preproc_anat_dir, [b '_coreg_gmseg.nii']);
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'flirt' ...
            ' -in ' PREPROC.anat_reference_file_gmseg ...
            ' -ref ' ref ...
            ' -applyxfm -init ' PREPROC.transform_T1_to_EPI{i} ...
            ' -interp nearestneighbour' ...
            ' -out ' PREPROC.coregistered_gmseg{i}]);
    end
    PREPROC.coregistered_wmseg_nuisance_ero{i, 1} = fullfile(PREPROC.preproc_anat_dir, [b '_coreg_wmseg_nuisance_ero.nii']);
    PREPROC.coregistered_csfseg_nuisance_ero{i, 1} = fullfile(PREPROC.preproc_anat_dir, [b '_coreg_csfseg_nuisance_ero.nii']);
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'flirt' ...
        ' -in ' PREPROC.anat_reference_file_wmseg_nuisance_erosion{end} ...
        ' -ref ' ref ...
        ' -applyxfm -init ' PREPROC.transform_T1_to_EPI{i} ...
        ' -interp nearestneighbour' ...
        ' -out ' PREPROC.coregistered_wmseg_nuisance_ero{i}]);
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'flirt' ...
        ' -in ' PREPROC.anat_reference_file_csfseg_nuisance_erosion{end} ...
        ' -ref ' ref ...
        ' -applyxfm -init ' PREPROC.transform_T1_to_EPI{i} ...
        ' -interp nearestneighbour' ...
        ' -out ' PREPROC.coregistered_csfseg_nuisance_ero{i}]);
    
    % Display coregistration result
    if do_check
        fprintf('Checking coregistration result...');
        system(['fsleyes -ad ' PREPROC.anat_reference_file ' ' PREPROC.coregistered_func_reference_files{i} ' &']);
    end
    
    close all;
    
end

% Take snapshot of coregistered reference images
fprintf('Taking snapshot of coregistered reference images.\n');
canlab_preproc_show_montage(PREPROC.coregistered_func_reference_files_masked, fullfile(PREPROC.qcdir, 'coreg_func_reference_masked.png'));
drawnow;

PREPROC = save_load_PREPROC(PREPROC.preproc_outputdir, 'save', PREPROC);
fprintf('\n\n\n');

end