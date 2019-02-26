function PREPROC = preproc_example(subjects)

% This function is an example of automatical preprocessing of CAPS2 dataset.
%
%
% :Usage:
% ::
%     PREPROC = preproc_example(subjects)
%
%
% :Input:
% ::
%   - subjects           subjects to be processed.
%                        (e.g., 1, 3, 18, [1:3], [1,2,18,20])
%
%
% :Output:
% ::
%     PREPROC
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

%% Get directory names and set basic environment

[rootdir, basedir, rscdir] = set_path_env('ncpu', 1);

%% Study directory setting and load data

study_imaging_dir = fullfile(basedir, 'data/CAPS2/Imaging/surface_analysis');

load(fullfile(basedir, 'projects/CAPS_project/data/CAPS2_dataset_171110.mat'));

%% Parameter setting

disdaq_n = 22; % Disdaq 22 images
tr = 460; % TR 460 msec
ciftify_basedir = fullfile(rscdir, 'fmri_toolboxes/CIFTIFY');
oasis_dir = fullfile(rscdir, 'OASIS_template');
n_thread = 1;
epi_enc_dir = 'ap';
ica_aroma_basedir = fullfile(rscdir, 'fmri_toolboxes/ICA-AROMA');
n_dim = 200;
fwhm = 5;
n_detrend = 1;
bandpass = [0 0.1];
wmcsf_method = 'mean';

%% Preprocessing!!!

for sj_num = subjects
    %% Subject information
    
    % e.g., 'CAPS2_002_AMK_170919' <- index for subject id: 7~9
    subject_code = sprintf('sub-caps%.3d', str2num(D.Subj_Level.id{sj_num}(7:9)));
    
    orderlist = D.Subj_Level.data(sj_num, 3:6);
    func_tasks = [];
    func_tasks{orderlist(1)} = 'CAPS';
    func_tasks{orderlist(2)} = 'QUIN';
    func_tasks{orderlist(3)} = 'ODOR';
    func_tasks{orderlist(4)} = 'REST';
    
    fprintf('study_imaging_dir: %s\nsubject_code: %s\n\n', study_imaging_dir, subject_code);
    disp(table(func_tasks));
    
    %% R-1. Make directories
    
    PREPROC = humanfmri_r1_make_directories(subject_code, study_imaging_dir, func_tasks);
    
    %% Copy DICOM images to new directory
    
    if ~exist('PREPROC', 'var')
        save_load_PREPROC(fullfile(study_imaging_dir, 'raw', subject_code), 'load');
    end
    
    print_header('Copy DICOM images to new directory: ', PREPROC.subject_code);
    
    imgdir = fullfile(basedir, 'data/CAPS2/Imaging/dicom_from_scanner', D.Subj_Level.id{sj_num});
    
    % source directories
    srcdir = [];
    srcdir.T1 = filenames(fullfile(imgdir, 'T1*'));
    srcdir.capsref = filenames(fullfile(imgdir, 'CAPS*SBREF*'));
    srcdir.caps = setdiff(filenames(fullfile(imgdir, 'CAPS*')), srcdir.capsref);
    srcdir.quinref = filenames(fullfile(imgdir, 'QUIN*SBREF*'));
    srcdir.quin = setdiff(filenames(fullfile(imgdir, 'QUIN*')), srcdir.quinref);
    srcdir.odorref = filenames(fullfile(imgdir, 'ODOR*SBREF*'));
    srcdir.odor = setdiff(filenames(fullfile(imgdir, 'ODOR*')), srcdir.odorref);
    srcdir.restref = filenames(fullfile(imgdir, 'REST*SBREF*'));
    srcdir.rest = setdiff(filenames(fullfile(imgdir, 'REST*')), srcdir.restref);
    srcdir.fmap_ap = filenames(fullfile(imgdir, 'DISTORTION*POLARITY_INVERT_TO_AP*'));
    srcdir.fmap_pa = setdiff(filenames(fullfile(imgdir, 'DISTORTION*')), srcdir.fmap_ap);
    
    % target directories
    dstndir = [];
    dstndir.T1 = PREPROC.dicom_anat_dir{1};
    dstndir.capsref = PREPROC.dicom_func_sbref_dir{contains(PREPROC.dicom_func_sbref_dir, 'task-CAPS')};
    dstndir.caps = PREPROC.dicom_func_bold_dir{contains(PREPROC.dicom_func_bold_dir, 'task-CAPS')};
    dstndir.quinref = PREPROC.dicom_func_sbref_dir{contains(PREPROC.dicom_func_sbref_dir, 'task-QUIN')};
    dstndir.quin = PREPROC.dicom_func_bold_dir{contains(PREPROC.dicom_func_bold_dir, 'task-QUIN')};
    dstndir.odorref = PREPROC.dicom_func_sbref_dir{contains(PREPROC.dicom_func_sbref_dir, 'task-ODOR')};
    dstndir.odor = PREPROC.dicom_func_bold_dir{contains(PREPROC.dicom_func_bold_dir, 'task-ODOR')};
    dstndir.restref = PREPROC.dicom_func_sbref_dir{contains(PREPROC.dicom_func_sbref_dir, 'task-REST')};
    dstndir.rest = PREPROC.dicom_func_bold_dir{contains(PREPROC.dicom_func_bold_dir, 'task-REST')};
    dstndir.fmap_ap = PREPROC.dicom_fmap_dir{contains(PREPROC.dicom_fmap_dir, 'dir-ap')};
    dstndir.fmap_pa = PREPROC.dicom_fmap_dir{contains(PREPROC.dicom_fmap_dir, 'dir-pa')};
    
    % check and copy
    fprintf('\n\nplease check source and destination directories.\n\n');
    fprintf('SOURCE directories:\n');
    disp(srcdir);
    fprintf('DESTINATION directories:\n');
    disp(dstndir);
    % input('to continue, press any key.\n\n');
    
    dirlist = fields(srcdir);
    for dir_i = 1:numel(dirlist)
        copy_src = char(getfield(srcdir, dirlist{dir_i}));
        copy_dstn = char(getfield(dstndir, dirlist{dir_i}));
        fprintf('\ncopying %s...\nSource: %s\nDestination: %s\n', dirlist{dir_i}, copy_src, copy_dstn);
        system(['cp -r ' copy_src ' ' copy_dstn]);
    end
    
    %% R-2. Dicom to NIFTI: structural
    
    PREPROC = humanfmri_r2_structural_dicom2nifti_bids(subject_code, study_imaging_dir);
    
    %% R-3. Dicom to nifti: functional
    
    PREPROC = humanfmri_r3_functional_dicom2nifti_bids(subject_code, study_imaging_dir, disdaq_n, 'tr', tr);
    
    %% R-4. Dicom to nifti: fmap
    
    PREPROC = humanfmri_r4_fieldmap_dicom2nifti_bids(subject_code, study_imaging_dir);
    
    %% S-1. Preproc directories
    
    PREPROC = humanfmri_s1_preproc_settings(subject_code, study_imaging_dir);
    
    %% S-2. Recon-all
    
    PREPROC = humanfmri_s2_surfrecon(subject_code, study_imaging_dir);
    
    %% S-3. CIFTIFY: ciftify_recon_all
    
    PREPROC = humanfmri_s3_ciftifysurf(subject_code, study_imaging_dir, ciftify_basedir);
    
    %% S-4. Anatomical segmentation and normalization
    
    % PREPROC = humanfmri_s4_anatomical_segmentation_normalization(subject_code, study_imaging_dir, oasis_dir, 'n_thread', n_thread);

    %% S-5. Estimate movement and correct misalignment between slices
    
    % PREPROC = humanfmri_s5_estimating_movement_slice_timing(subject_code, study_imaging_dir);
    
    %% S-6. Motion correction
    
    PREPROC = humanfmri_s6_motion_correction(subject_code, study_imaging_dir);
    
    %% S-7. Distortion correction
    
    PREPROC = humanfmri_s7_distortion_correction(subject_code, study_imaging_dir, 'epi_enc_dir', epi_enc_dir);
    
    %% S-8. Coregistration
    
    PREPROC = humanfmri_s8_coregistration(subject_code, study_imaging_dir, 'no_check_reg');
    
    %% S-9. ICA-AROMA
    
    PREPROC = humanfmri_s9_ICA_AROMA(subject_code, study_imaging_dir, ica_aroma_basedir, 'n_dim', n_dim, 'n_thread', n_thread);
    
    %% S-10. Denoising
    
    PREPROC = save_load_PREPROC(fullfile(study_imaging_dir, 'preprocessed', subject_code), 'load'); % load PREPROC
    rating_cov = [];
    for run_i = 1:4
        disdaq_correction = PREPROC.func_disdaq_number(run_i) * PREPROC.func_TR(run_i) / 1000 - 10; % sec
        
        wh_change = logical(D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Changecolor_Stim'])});
        wh_abruptchange = find(diff([0;wh_change])~=0);
        wh_change_starttime = D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Timefromstart'])}(wh_abruptchange(1:2:end));
        wh_change_positive = wh_change_starttime - disdaq_correction > 0;
        wh_change_starttime = wh_change_starttime(wh_change_positive);
        wh_change_endtime = D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Timefromstart'])}(wh_abruptchange(2:2:end));
        wh_change_endtime = wh_change_endtime(wh_change_positive);
        wh_change_duration = wh_change_endtime - wh_change_starttime;
        
        wh_response = logical(D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Changecolor_Response'])});
        wh_abruptresponse = find(diff([0;wh_response])~=0);
        wh_response_starttime = D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Timefromstart'])}(wh_abruptresponse(1:2:end));
        wh_response_positive = wh_response_starttime - disdaq_correction > 0;
        wh_response_starttime = wh_response_starttime(wh_response_positive);
        wh_response_endtime = D.Continuous.data{sj_num}{strcmp(D.Continuous.names, ['Run', num2str(run_i), '_Timefromstart'])}(wh_abruptresponse(2:2:end));
        wh_response_endtime = wh_response_endtime(wh_response_positive);
        wh_response_duration = wh_response_endtime - wh_response_starttime;
        
        rating_cov{run_i} = [onsets2fmridesign([wh_change_starttime - disdaq_correction, wh_change_duration], PREPROC.func_TR(run_i)/1000, PREPROC.func_nvol(run_i)*PREPROC.func_TR(run_i)/1000, spm_hrf(1)), ...
            onsets2fmridesign([wh_response_starttime - disdaq_correction, wh_response_duration], PREPROC.func_TR(run_i)/1000, PREPROC.func_nvol(run_i)*PREPROC.func_TR(run_i)/1000, spm_hrf(1))];
        rating_cov{run_i} = rating_cov{run_i}(:,[1,3]); % not intercept
    end
    
    PREPROC = humanfmri_s10_denoising(subject_code, study_imaging_dir, 'detrend', n_detrend, 'bandpass', bandpass, 'wmcsf', wmcsf_method, 'custom_nuisance', rating_cov, 'addmean');
    
    %% S-11. Normalization
    
    PREPROC = humanfmri_s11_normalization(subject_code, study_imaging_dir, 'n_thread', n_thread);
    
    %% S-12. Smoothing
    
    PREPROC = humanfmri_s12_smoothing(subject_code, study_imaging_dir, 'fwhm', fwhm);
    
    %% S-13. CIFTIFY: ciftify_subject_fmri
    
    PREPROC = humanfmri_s13_ciftifyfmri(subject_code, study_imaging_dir, ciftify_basedir, 'fwhm', fwhm);
    
end

end
