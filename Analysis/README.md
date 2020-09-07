# Zcalbase-gal/Analysis

1. Overview
2. Execution 
    1. Requirements
    2. Running Analysis 
3. Changelog
4. Authors

# Overview

Python 3.7.7 codes analysis of a subset of the [DEEP2 Redshift Survey](http://deep.ps.uci.edu).
Code relies on [Metallicity Stack Commons](https://github.com/astrochun/Metallicity_Stack_Commons)
created by Chun Ly, Reagen Leimbach, and Caroline McCormick. 

Previous studies of galaxy evolution have been proven incorrect due to the fact that the study did not take into 
account the difference in the gas properties between high- and low-z galaxies (redshifted). This analysis presents
a calibration between the metallicity and strong-line diagnostics, R23 and O32, to use to correct for the change 
in gas properties as redshift increases. To create this calibration, the OIII$\lambda$[4363] is measured to calculate
the electron temperature of the galaxy, which is used to calculate the metallicity. 

OIII$\lambda$[4363] is a forbidden transition, which produces a very weak emission line in the spectra. In order to 
detect OIII$\lambda$[4363], individual spectra are binned and stacked to get a composite spectrum measurement. 

This analysis is also compared against recently published analyses from 
Jiang, T., Malhotra, S., Rhoads, J. E., et al. 2019, ApJ, 872, 145  and 
Bian, F., Kewley, L. J., & Dopita, M. A. 2018, ApJ, 859, 175. 

Links to material presenting  to 


# Execution 

### Requirements 
Your will need the following to have a working copy of this software.

- Python(3.7)
- [Metallicity Stack Commons](https://github.com/astrochun/Metallicity_Stack_Commons)
- numpy 
- matplotlib
- astropy
- scipy

### Running Analysis 
The analysis is run by executing the run_grid_R23_O32_analysis() function in Analysis/DEEP2_R23_O32/general.py. 
By importing analysis/deep2_r23_o32/general.py in a python environment, the codes needed for the analysis 
(see list below) are imported.


Once run_grid_R23_O32_analysis() has run, other codes not listed below can be used to create various plots
in addition what is created during the run. 


