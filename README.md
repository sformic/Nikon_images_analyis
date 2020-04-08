# Nikon images analysis 
### Sara Formichetti

## Aim

Understanding whether there is a difference in the nuclear localization of OGT protein (and consequently of O-GlcNAc modification) between wild-type MEFs and MEFs derived from male mouse embryos whose endogenous Ogt gene bears a mutation in a putative Nuclear Localization Signal (NLS) according to [Seo et al. Sci Rep 2016](https://www.nature.com/articles/srep34614).

## Experimental design

Immunofluorescence with antibodies against OGT, O-GlcNAc and OGA.

Cell model: MEFs

Samples:

* wild-type MEFs (1 clone), one well per antibody
* NLSm-OGT MEFs (i.e. MEFs bearing a variant of endogenous OGT with a mutation in a putative Nuclear Localization Signal of OGT) (1 clone), one well per antibody

Imaging:

Microscope: Nikon A1 confocal microscope in trial at EMBL Rome in December 2019.
Strategy: Per sample per antibody: big picture (5x5 60x with Galvano laser) --> z-stack of all picture (from 8 to 11 depending on the picture) --> output is nd2 file to open with fiji and analyze a posteriori (see below) 

## Image analysis

General idea: finding "central" nuclear plane for 30 cells --> quantifying OGT/O-GlcNAc fluorescence at that plane --> checking statistically significant difference between wt and mutant MEFs

### Practical steps

0. I could not load images to github neither from brutus to scrap because too big so I copied one image (5 ab177941.nd2 i.e. antiOGT antibody, clone 5) on the Desktop in the boulard account of Imaris workstation. pw for boulard user is microboulard
1.  [Macro file](./scripts/Macro.ijm) which, for each .nd2 image:

* opens the DAPi channel and the other channels separately
* bluring, binarization, Distance Transform Watershed 3D --> label mask image 3D
* "stack to images" for both intensity original stack and label mask image 3D 
* measure of mean, median etc for each plane using corresponding plane's label mask image

3. [R script](./scripts/nuclei_param_change_across_planes.R) which I did in order to understand what should be the best parameter to use to choose the central nuclear plane. Here I stopped because some plots made me noticing I have a problem: I will explain you better in the VC.
So far I arrived here. After that I was thinking to take, for each label i.e. nucleus, only the "central" plane, then measure OGT/O-GlcNAc intensity there. This part should be more straightforward. 

## My issues and doubts until this point:

1. With Distance Transform Watershed I did not manage to find the parameters to separate attached nuclei. DIfferntly from what I found in the internet as examples for nuclei segmentation, in my case I have little nuclei distant to each other with some cases of doublets. I am afraid this is difficult for watershed but maybe I am wrong.
2. Why is Intensity measurements command measuring the VOLUME if I gave it 2D images to compare??
3. Ideally, The best in my opinion would be to measure the mean DAPI intensity in each plane but using the 3D nuclear segmentation as label mask. However, I did not find a practical way to do so. 
