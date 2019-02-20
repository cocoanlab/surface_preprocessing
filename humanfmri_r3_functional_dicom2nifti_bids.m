function PREPROC = humanfmri_r3_functional_dicom2nifti_bids(subject_code, study_imaging_dir, disdaq_n, varargin)

% This function convert the functional dicom files into nifti and json files.
%
%
% :Usage:
% ::
%     PREPROC = humanfmri_r3_functional_dicom2nifti_bids(subject_code, study_imaging_dir, disdaq_n, varargin)
%
%
% :Input:
% ::
%   - subject_code       the subject code (e.g., 'sub-caps001').
%   - study_imaging_dir  the directory information for the study imaging data 
%                        (e.g., '/Volumes/habenula/hbmnas/data/CAPS2/Imaging').
%   - disdaq_n           the number of images you want to discard (to allow for
%                        image intensity stablization) for functional runs; it
%                        can be one value or multiple values corresponding
%                        to each run.
%                        (e.g., disdaq_n = 22
%                               -> discard 22 images for every run,
%                               disdaq_n = [22 20 21]
%                               -> discard 22 images for run 1,
%                               -> discard 20 images for run 2,
%                               -> discard 21 images for run 3)
%            
%
% :Optional Input:
% ::
%   - run_num            runs to include. ex) [1 2 4 5], ...
%   - tr                 repetition time (in milliseconds).
%                        (e.g., 460
%                               -> TR = 460 msec for every run,
%                               [460, 850, 780]
%                               -> TR = 460 msec for run 1,
%                               -> TR = 850 msec for run 2,
%                               -> TR = 780 msec for run 3)
%                        if not specified, it is obtained by reading dicom
%                        header file.
%
%
% :Output:
% ::
%     PREPROC.raw_func_dir
%     PREPROC.func_bold_files
%     PREPROC.func_bold_json_files
%     PREPROC.func_bold_dicomheader_files
%     PREPROC.func_sbref_files
%     PREPROC.func_sbref_json_files
%     PREPROC.func_sbref_dicomheader_files
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
tr = NaN;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'run_num'}
                do_select_run = true;
                run_num = varargin{i+1};
            case {'tr'}
                tr = varargin{i+1};
        end
    end
end


subject_dir = fullfile(study_imaging_dir, 'raw', subject_code);
PREPROC = save_load_PREPROC(subject_dir, 'load'); % load PREPROC

print_header('Convert FUNCTIONAL DICOM images to NIFTI format', PREPROC.subject_code);

dcmheaders_dir = fullfile(PREPROC.study_imaging_dir, 'dcmheaders', PREPROC.subject_code);
system(['mkdir -p ' dcmheaders_dir]);

PREPROC.raw_func_dir = fullfile(PREPROC.subject_dir, 'func');
system(['mkdir -p ' PREPROC.raw_func_dir]);

if numel(disdaq_n) == 1
    disdaq_n = repmat(disdaq_n, numel(PREPROC.dicom_func_bold_dir), 1);
end
disdaq_n = disdaq_n(:);
if numel(tr) == 1
    tr = repmat(tr, numel(PREPROC.dicom_func_bold_dir), 1);
end
tr = tr(:);

for i = 1:numel(PREPROC.dicom_func_bold_dir)
    
    if ~do_select_run || ismember(i, run_num)
        
        [~, functype] = fileparts(PREPROC.dicom_func_bold_dir{i});
        [regexp_start, regexp_end] = regexp(functype, 'task-\w*_run-');
        taskname = functype(regexp_start+length('task-'):regexp_end-length('_run-'));
        fprintf('\n\nFunctional image type: %s\n\n', taskname);
        
        dicom_imgs = search_files(fullfile(PREPROC.dicom_func_bold_dir{i}, '*.IMA'), 3);
        [~, temp_dicom_dir] = system('mktemp -d');
        temp_dicom_dir = strtrim(temp_dicom_dir);
        
        dicm2nii(dicom_imgs, temp_dicom_dir, 4, 'save_json', 'taskname', taskname);
        
        load(fullfile(temp_dicom_dir, 'dcmHeaders.mat'));
        dicm2nii_name = char(fieldnames(h));
        h = getfield(h, dicm2nii_name);
        PREPROC.func_bold_dicomheader_files{i, 1} = fullfile(dcmheaders_dir, [PREPROC.subject_code '_' functype '_dcmheaders.mat']);
        save(PREPROC.func_bold_dicomheader_files{i}, 'h');
        system(['rm ' fullfile(temp_dicom_dir, 'dcmHeaders.mat')]);
        
        PREPROC.func_disdaq_number(i, 1) = disdaq_n(i);
        if PREPROC.func_disdaq_number(i) > 0
            fprintf('The number of images to be discarded: %d\n', PREPROC.func_disdaq_number(i));
            system(['cd ' temp_dicom_dir ';' ...
                'dicm2nii_name=($ls ' [dicm2nii_name '*.nii'] ');' ...
                'for i in {0..' num2str(PREPROC.func_disdaq_number(i)-1) '}; do ' ...
                'rm ${dicm2nii_name[$i]}' ...
                '; done']);
        end
        
        if isnan(tr(i))
            tr(i) = h.RepetitionTime; % msec
        end
        PREPROC.func_TR(i, 1) = tr(i);
        fprintf('TR: %d msec\n', PREPROC.func_TR(i));
        
        PREPROC.func_bold_files{i, 1} = fullfile(PREPROC.raw_func_dir, [PREPROC.subject_code '_' functype '.nii']);
        fprintf('Merging 3d images to 4d images...\n')
        system(['export FSLOUTPUTTYPE=NIFTI;' ...
            ...
            'cd ' temp_dicom_dir ';' ...
            'fslmerge' ...
            ' -tr' ...
            ' ' PREPROC.func_bold_files{i} ...
            ' ' [dicm2nii_name '*.nii'] ...
            ' ' num2str(PREPROC.func_TR(i) / 1000)]); % sec
        
        [~, sliceinfo] = system(['echo $(fslval ' PREPROC.func_bold_files{i} ' dim3)']);
        PREPROC.func_nslices(i, 1) = str2num(sliceinfo);
        fprintf('The number of slices: %d\n', PREPROC.func_nslices(i));
        [~, volinfo] = system(['echo $(fslval ' PREPROC.func_bold_files{i} ' dim4)']);
        PREPROC.func_nvol(i, 1) = str2num(volinfo);
        fprintf('The number of volumes: %d\n', PREPROC.func_nvol(i));
        
        PREPROC.func_bold_json_files{i, 1} = fullfile(PREPROC.raw_func_dir, [PREPROC.subject_code '_' functype '.json']);
        system(['mv ' fullfile(temp_dicom_dir, [dicm2nii_name '.json']) ' ' PREPROC.func_bold_json_files{i}]);
        
        system(['rm -r ' temp_dicom_dir]);
        
        if ~isempty(PREPROC.dicom_func_sbref_dir{i})
            
            [~, functype] = fileparts(PREPROC.dicom_func_sbref_dir{i});
            
            dicom_imgs = search_files(fullfile(PREPROC.dicom_func_sbref_dir{i}, '*.IMA'), 3);
            
            dicm2nii(dicom_imgs, PREPROC.raw_func_dir, 4, 'save_json', 'taskname', taskname);
            
            load(fullfile(PREPROC.raw_func_dir, 'dcmHeaders.mat'));
            dicm2nii_name = char(fieldnames(h));
            h = getfield(h, dicm2nii_name);
            PREPROC.func_sbref_dicomheader_files{i, 1} = fullfile(dcmheaders_dir, [PREPROC.subject_code '_' functype '_dcmheaders.mat']);
            save(PREPROC.func_sbref_dicomheader_files{i}, 'h');
            system(['rm ' fullfile(PREPROC.raw_func_dir, 'dcmHeaders.mat')]);
    
            PREPROC.func_sbref_files{i, 1} = fullfile(PREPROC.raw_func_dir, [PREPROC.subject_code '_' functype '.nii']);
            system(['mv ' fullfile(PREPROC.raw_func_dir, [dicm2nii_name '.nii']) ' ' PREPROC.func_sbref_files{i}]);
            
            PREPROC.func_sbref_json_files{i, 1} = fullfile(PREPROC.raw_func_dir, [PREPROC.subject_code '_' functype '.json']);
            system(['mv ' fullfile(PREPROC.raw_func_dir, [dicm2nii_name '.json']) ' ' PREPROC.func_sbref_json_files{i}]);
            
        end
        
    end

end

save_load_PREPROC(PREPROC.subject_dir, 'save', PREPROC); % save PREPROC
fprintf('\n\n\n');

end