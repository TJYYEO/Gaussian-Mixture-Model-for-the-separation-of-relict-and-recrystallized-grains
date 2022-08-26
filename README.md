# Gaussian-Mixture-Model-for-the-separation-of-relict-and-recrystallized-grains
[![DOI](https://zenodo.org/badge/529074859.svg)](https://zenodo.org/badge/latestdoi/529074859)

Matlab script for processing electron backscatter diffraction (EBSD) data. 

We advise for the calculations to be performed in MATLAB R2021b and MTEX 5.7.0 or any available later versions.

‘EBSD_Grains_Construction.mlx’ mainly consists of procedures for noise removal and restrictions applied when constructing the grains. Noise includes wild spikes and systematic mis-indexed pixels. Dauphine twins are merged with their parent grains.

‘RelictRex_Separator.mlx’ and ‘GaussianMM.m’ are used for the evaluation of relict and recrystallized grains.

Files (15.9 kB)
