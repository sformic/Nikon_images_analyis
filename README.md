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

* Microscope: Nikon A1 confocal microscope in trial at EMBL Rome in December 2019.
* Strategy: Per sample per antibody: big picture (5x5 60x with Galvano laser) --> z-stack of all picture (from 8 to 11 depending on the picture) --> output is nd2 file to open with fiji and analyze a posteriori (see below) 

## Image analysis

General idea: finding "central" nuclear plane for 30 cells --> quantifying OGT/O-GlcNAc fluorescence at that plane --> checking statistically significant difference between wt and mutant MEFs.
Controlling for generalized (i.e. both in cytosol and nucleus) staining difference between wt and mutant sample by comparing intensities in the cytosolic areas. 

### Practical steps

0. [Rmd that uses manual measures](./scripts/groundtruth_measures.Rmd) inside manually drawn areas to verify that the strategy of measuring ch1 at the plane with highest DAPI intensity makes sense. Download the relative html to see outcomes and plots.

1.  [Macro file](./scripts/Macro_2.ijm), which should be commented enough to understand all steps and rationales. This Macro is very different from [the initial one](./scripts/Macro.ijm) and it is based on last suggestions by Alvaro.

2. [R script](./scripts/ch_at_best_nuclear_plane.R) that uses all tables output by the Macro_2 algorithm and test statistical significance of difference of ch1 nuclear intensity between wt and mutant MEFs for all antibodies in folder smb://brutus.embl.it/boulard/Sara/Julia/Data/confocal/191216 Sara-Julia Nikon A1.

The final plots resulting from the analysis are in [this folder](./scripts/images).