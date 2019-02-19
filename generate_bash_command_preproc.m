%%
[rootdir, basedir, rscdir] = set_path_env('ncpu', 1);
load(fullfile(basedir, 'projects/CAPS_project/data/CAPS2_dataset_171110.mat'));

%%
sj_num = 8;
run_preproc_command = ['matlab -nodesktop -nosplash -nodisplay -r "addpath(genpath(''' fullfile(rscdir, 'matlab_toolboxes/surface_preprocessing') '''));' ...
    'CAPS2_preproc([' num2str(sj_num) ']); quit" >> ~/' sprintf('sub-caps%.3d_log.txt', sj_num) ' 2>&1 < /dev/null &'];
clipboard('copy', run_preproc_command);

