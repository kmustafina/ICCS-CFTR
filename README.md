# ICCS-CFTR
This repo contains Wiseman Lab's MATLAB implementation of ICCS for the paper

## Documentation
The "functions" folder contains custom MATLAB functions needed to perform the analysis.

- autocrop.m: automatic cropping of the correlation function
- corrfunc.m: calculation of autocorrelation functions
- corrfunc_cross.m: calculation of cross-correlation functions
- gauss2d.m, gauss2dwxy.m, gaussfit.m: Gaussian fitting of the correlation functions
- wnCorr.m: background correction for images

Analysis can be tested by running the runfile_example.m on the bead data in the "data" folder.

## Requirements
MATLAB 2020b

## Acknowledgements
The majority of code was written by David Kolin, a previous Wiseman Lab member.
