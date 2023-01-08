# Gaussian-Mixture-Model-for-the-separation-of-relict-and-recrystallized-grains v1.0
[![DOI](https://zenodo.org/badge/529074859.svg)](https://zenodo.org/badge/latestdoi/529074859)

# Gaussian-Mixture-Model-for-the-separation-of-relict-and-recrystallized-grains v1.1a
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7513423.svg)](https://doi.org/10.5281/zenodo.7513423)

Matlab script for processing electron backscatter diffraction (EBSD) data. 

We advise for the calculations to be performed in MATLAB R2021b and MTEX 5.7.0 or any available later versions.

‘EBSD_Grains_Construction.mlx’ mainly consists of procedures for noise removal and restrictions applied when constructing the grains. Noise includes wild spikes and systematic mis-indexed pixels. Dauphine twins are merged with their parent grains.

‘RelictRex_Separator.mlx’ and ‘GaussianMM.m’ are used for the evaluation of relict and recrystallized grains.

Version 1.1 has been updated to optimised the number of iterations necessary for the Monte Carlo simulation.

Files (15.9 kB)
