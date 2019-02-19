function PREPROC = humanfmri_r2_anatomical_dicom2nifti_bids(subject_code, study_imaging_dir)

% This function convert the anatomical dicom files into nifti and json files.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_r2_anatomical_dicom2nifti_bids(subject_code, study_imaging_dir)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%
%
% :Output:
% ::
%     PREPROC.raw_anat_dir
%     PREPROC.anat_files
%     PREPROC.anat_json_files
%     PREPROC.anat_dicomheader_files
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
subject_dir = fullfile(study_imaging_dir, 'raw', subject_code);
PREPROC = save_load_PREPROC(subject_dir, 'load'); % load PREPROC

print_header('Convert ANATOMICAL DICOM images to NIFTI format', PREPROC.subject_code);

dcmheaders_dir = fullfile(PREPROC.study_imaging_dir, 'dcmheaders', PREPROC.subject_code);
system(['mkdir -p ' dcmheaders_dir]);

PREPROC.raw_anat_dir = fullfile(PREPROC.subject_dir, 'anat');
system(['mkdir -p ' PREPROC.raw_anat_dir]);

for i = 1:numel(PREPROC.dicom_anat_dir)
    
    [~, anattype] = fileparts(PREPROC.dicom_anat_dir{i});
    fprintf('\n\nAnatomical image type: %s\n\n', anattype);
    
    dicom_imgs = search_files(fullfile(PREPROC.dicom_anat_dir{i}, '*.IMA'), 3);
    
    dicm2nii(dicom_imgs, PREPROC.raw_anat_dir, 4, 'save_json');
    
    load(fullfile(PREPROC.raw_anat_dir, 'dcmHeaders.mat'));
    dicm2nii_name = char(fieldnames(h));
    h = getfield(h, dicm2nii_name);
    PREPROC.anat_dicomheader_files{i, 1} = fullfile(dcmheaders_dir, [PREPROC.subject_code '_' anattype '_dcmheaders.mat']);
    save(PREPROC.anat_dicomheader_files{i}, 'h');
    system(['rm ' fullfile(PREPROC.raw_anat_dir, 'dcmHeaders.mat')]);
    
    PREPROC.anat_files{i, 1} = fullfile(PREPROC.raw_anat_dir, [PREPROC.subject_code '_' anattype '.nii']);
    system(['mv ' fullfile(PREPROC.raw_anat_dir, [dicm2nii_name '.nii']) ' ' PREPROC.anat_files{i}]);
    
    PREPROC.anat_json_files{i, 1} = fullfile(PREPROC.raw_anat_dir, [PREPROC.subject_code '_' anattype '.json']);
    system(['mv ' fullfile(PREPROC.raw_anat_dir, [dicm2nii_name '.json']) ' ' PREPROC.anat_json_files{i}]);
    
end

save_load_PREPROC(PREPROC.subject_dir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end