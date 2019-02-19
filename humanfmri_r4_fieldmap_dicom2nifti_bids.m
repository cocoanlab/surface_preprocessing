function PREPROC = humanfmri_r4_fieldmap_dicom2nifti_bids(subject_code, study_imaging_dir)

% This function convert the fieldmap dicom files into nifti and json files.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_r4_fieldmap_dicom2nifti_bids(subject_code, study_imaging_dir)
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
%     PREPROC.raw_fmap_dir
%     PREPROC.fmap_files
%     PREPROC.fmap_json_files
%     PREPROC.fmap_dicomheader_files
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

print_header('Convert FIELDMAP DICOM images to NIFTI format', PREPROC.subject_code);

dcmheaders_dir = fullfile(PREPROC.study_imaging_dir, 'dcmheaders', PREPROC.subject_code);
system(['mkdir -p ' dcmheaders_dir]);

PREPROC.raw_fmap_dir = fullfile(PREPROC.subject_dir, 'fmap');
system(['mkdir -p ' PREPROC.raw_fmap_dir]);

for i = 1:numel(PREPROC.dicom_fmap_dir)
    
    [~, fmaptype] = fileparts(PREPROC.dicom_fmap_dir{i});
    fprintf('\n\nFieldmap image type: %s\n\n', fmaptype);
    
    dicom_imgs = search_files(fullfile(PREPROC.dicom_fmap_dir{i}, '*.IMA'), 3);
    [~, temp_dicom_dir] = system('mktemp -d');
    temp_dicom_dir = strtrim(temp_dicom_dir);
    
    dicm2nii(dicom_imgs, temp_dicom_dir, 4, 'save_json');
        
    load(fullfile(temp_dicom_dir, 'dcmHeaders.mat'));
    dicm2nii_name = char(fieldnames(h));
    h = getfield(h, dicm2nii_name);
    PREPROC.fmap_dicomheader_files{i, 1} = fullfile(dcmheaders_dir, [PREPROC.subject_code '_' fmaptype '_dcmheaders.mat']);
    save(PREPROC.fmap_dicomheader_files{i}, 'h');
    system(['rm ' fullfile(temp_dicom_dir, 'dcmHeaders.mat')]);
    
    PREPROC.fmap_files{i, 1} = fullfile(PREPROC.raw_fmap_dir, [PREPROC.subject_code '_' fmaptype '.nii']);
    disp('Merging 3d images to 4d images...')
    system(['export FSLOUTPUTTYPE=NIFTI;' ...
        ...
        'cd ' temp_dicom_dir ';' ...
        'fslmerge' ...
        ' -t' ...
        ' ' PREPROC.fmap_files{i} ...
        ' ' [dicm2nii_name '*.nii']]);
    
    PREPROC.fmap_json_files{i, 1} = fullfile(PREPROC.raw_fmap_dir, [PREPROC.subject_code '_' fmaptype '.json']);
    system(['mv ' fullfile(temp_dicom_dir, [dicm2nii_name '.json']) ' ' PREPROC.fmap_json_files{i}]);
    
    system(['rm -r ' temp_dicom_dir]);

end

save_load_PREPROC(PREPROC.subject_dir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end