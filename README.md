# Surface Preprocessing Pipeline for human fMRI (cocoanlab)

This repository includes a set of matlab functions for fMRI data preprocessing. Particularly, this pipeline was made for surface-based preprocessing, but it also can be used for volume-based preprocessing.

Overall scheme of this pipeline is as follows:

A) Make directories for preprocessed data (r1)

B) Convert DICOM images to NIFTI format (using dicm2nii.m, which was adapted from https://github.com/xiangruili/dicm2nii) (r2-r4)

C) Basic environmental setup for preprocessing (s1)

D) Preprocessing structural data

  For surface-based preprocessing:
  
  1) Use recon-all (Freesurfer) to correct bias-field, extract brain tissue, reconstruct structural images to cortical surface, and do anatomical segmentation (s2)
  2) Use ciftify_recon_all (CIFTIFY, Dickie et al., 2019) to normalize structural images to MNI space, and resample native surface images onto Conte69 164k and 32k CIFTI surface (Van Essen et al., 2012) using MNI normalization parameters (s3)
    
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
  6) Use ICA-AROMA (Pruim et al., 2015) to remove motion-related signals (s9)
  7) Use 3dTproject (AFNI) to remove nuisance signals (s10)
  8) Use applywarp (FSL, for surface-based preprocessing) or antsApplyTransforms (ANTs, for volume-based preprocessing) to normalize functional images to MNI space
  9) For volume-based preprocessing, use susan (FSL) to spatially smooth functional images
  10) For surface-based preprocessing, use ciftify_subject_fmri (CIFTIFY, Dickie et al., 2019) to transform functional images to Conte69 32k CIFTI surface, and spatially smooth (both surface-based and volume-based smoothing) functional images
  
  Note: 8) can precede 6) and 7). depending on your purpose.


Dependency:
  CanlabCore: [https://github.com/canlab/CanlabCore](https://github.com/canlab/CanlabCore)
  dicm2nii: ([link](https://www.mathworks.com/matlabcentral/fileexchange/42997-dicom-to-nifti-converter--nifti-tool-and-viewer)): We modified the original toolbox a little bit to make the output data fully BIDS-compatible. For this reason, please use the dicm2nii toolbox in our repository (in /external), instead of the original one. 
  AFNI: [https://afni.nimh.nih.gov](https://afni.nimh.nih.gov)
  FSL: [https://fsl.fmrib.ox.ac.uk](https://fsl.fmrib.ox.ac.uk)
  ANTs: [https://github.com/ANTsX/ANTs](https://github.com/ANTsX/ANTs)
  Freesurfer: [https://surfer.nmr.mgh.harvard.edu](https://surfer.nmr.mgh.harvard.edu)
  Connectome Workbench: [https://www.humanconnectome.org/software/connectome-workbench](https://www.humanconnectome.org/software/connectome-workbench)
  ICA-AROMA: [https://github.com/rhr-pruim/ICA-AROMA](https://github.com/rhr-pruim/ICA-AROMA)
  
  
2019. 07. 20
Jae-Joong Lee
