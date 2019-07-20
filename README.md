# Surface Preprocessing Pipeline for human fMRI (cocoanlab)

This repository includes a set of matlab functions for fMRI data preprocessing. Particularly, this pipeline was made for surface-based preprocessing, but it also can be used for volume-based preprocessing.

Overall scheme of this pipeline is as follows:

A) Make directories for preprocessed data (r1)

B) Convert DICOM images to NIFTI format (using dicm2nii.m, which was adapted from https://github.com/xiangruili/dicm2nii) (r2-r4)

C) Basic environmental setup for preprocessing (s1)

D) Preprocessing structural data

  For surface-based preprocessing:
  
  1) Use recon-all (Freesurfer) to correct bias-field, extract brain tissue, reconstruct structural images to cortical surface, and do anatomical segmentation (s2)
  2) Use ciftify_recon_all (CIFTIFY, Dickie et al., 2019) to normalize structural images to MNI space, and resample native surface images onto Conte69 164k and 32k surface (Van Essen et al., 2012) using MNI normalization parameters (s3)
    
  For Volume-based preprocessing:
    
  1) Use antsBrainExtraction (ANTs) to correct bias-field and extract brain tissue (s4)
  2) Use antsRegistrationSyN (ANTs) to normalize structural images to MNI space (s4)
  3) Use FAST (FSL) to do anatomical segmentation (s4)
  
E) Preprocessing functional data

  1) Use 3dTshift (AFNI) to do slice-timing correction (s5, if needed)
  2) Use 3dvolreg (AFNI) to do motion-correction (s6)
  3) Use topup/applytopup (FSL) to do distortion-correction (s7)
  4) Use flirt (FSL) with BBR cost function to co-register functional images to structural images (s8)
  5) For surface-based preprocessing, additionally use bbregster (Freesurfer) to refine co-registration (s8)
  6) 


uses a lot of useful softwares (AFNI, FSL, SPM, ANTs, freesurfer, Connectome Workbench, ICA-AROMA, and CIFTIFY)
to preprocess human fMRI datasets.

This repository includes a set of matlab functions for fMRI data preprocessing. Currently, this includes tools for a) dicom to nifti in the BIDS (Brain Imaging Data Structure) format, using dicm2nii (which has been a little bit modified, and therefore it is included in this toolbox), b) distortion correction using fsl's topup, c) outlier detection and data quality check using some useful Canlab tools (github), d) slice-timing, realignment, T1 (default) or EPI normalization, and smoothing using SPM12 tools, and e) a recent tool for motion correction, ICA-AROMA

Jae-Joong Lee
