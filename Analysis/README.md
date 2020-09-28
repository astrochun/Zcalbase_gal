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
in gas properties as redshift increases. To create this calibration, the OIIIλ[4363] is measured to calculate
the electron temperature of the galaxy, which is used to calculate the metallicity. 

OIIIλ[4363] is a forbidden transition, which produces a very weak emission line in the spectra. In order to 
detect OIIIλ[4363], individual spectra are binned and stacked to get a composite spectrum measurement. 

This analysis is also compared against recently published analyses from 
Jiang, T., Malhotra, S., Rhoads, J. E., et al. 2019, ApJ, 872, 145  and 
Bian, F., Kewley, L. J., & Dopita, M. A. 2018, ApJ, 859, 175. 

Links to material presenting this work: 
1. Link to figure created in the stacking code: 
https://arizona.figshare.com/articles/Stacking_of_Galaxy_Spectra/12360626

2. Link to Honors Thesis

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



    from Zcalbase_gal.analysis.deep2_r23_o32 import general


By importing analysis/deep2_r23_o32/general.py in a python environment, the following codes from Zcalbase_gal 
and other sub-folders are imported. 
- stackboth_mastergrid
- zoom_and_gauss_general
- hstack_tables
- r_temp_calcul
- calibration_plots
- binning/n_bins_grid_analysis
- binning/fixed_grid_analysis
- binning/single_grid_o32 
- binning/single_grid_r23
- plotting/more_plots
- plotting/line_ratio_plotting
- plotting/te_metal_plots

The run function requires the following variables. 

- dataset -> keyword used to define binning method  options: Grid, O32_Grid, R23_Grid, n_Bins
- y_correction -> determines if the smoothed (movingaverage_box1D) version of y is used in zoom_and_gauss_general.py
- n_split -> determined how many times the R23 bins are split when using manual binning
- adaptive -> determines if the R23 bins have equal or different number of spectra in them in binning method
- dustatten -> determines if dust attenuation corrections are applied
- mask -> determines if the night sky mask is used in Stackingboth_MasterGrid.py

Calling the run function


    general.run_grid_r23_o32_analysis(dataset, y_correction, n_split, 
    adaptive, dustatten, mask)
    

Steps taking throughout run function: 

a. 
