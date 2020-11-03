# Executing Zcalbase_gal on DEEP2 data

1. [Overview](#overview)
2. [Execution](#execution)
    1. [Requirements](#requirements)
    2. [Running Grid Analysis](#running-grid-analysis)
    3. [Running Voronoi Analysis](#running-voronoi-analysis)

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
[Jiang, T., Malhotra, S., Rhoads, J. E., et al. 2019, ApJ, 872, 145](https://arxiv.org/abs/1811.05796)  
and  [Bian, F., Kewley, L. J., & Dopita, M. A. 2018, ApJ, 859, 175](https://iopscience.iop.org/article/10.3847/1538-4357/aabd74/meta). 

Links to material presenting this work: 
1. [UA Space Grant Link](https://arizona.figshare.com/articles/Stacking_of_Galaxy_Spectra/12360626) 
to figure created in the stacking code 


2. [Link (TBD)] to Honors Thesis

# Execution 

### Requirements 
Your will need the following to have a working copy of this software.

- Python(3.7)
- [Metallicity Stack Commons](https://github.com/astrochun/Metallicity_Stack_Commons)
- numpy 
- matplotlib
- astropy
- scipy
- [chun_codes](https://github.com/astrochun/chun_codes)

### Running Grid Analysis 
The analysis of the binning methods is run by executing the `run_grid_R23_O32_analysis`
function in [Analysis/DEEP2_R23_O32/general.py](deep2_r23_o32/general.py). 
```python
    from Zcalbase_gal.analysis.deep2_r23_o32 import general
```
The `run` function requires the following variables.

- `dataset`: keyword used to define binning method  options: Grid, O32_Grid, R23_Grid, n_Bins
- `y_correction`: determines if the smoothed (movingaverage_box1D) version of y is used in zoom_and_gauss_general.py
- `n_split`: determined how many times the R23 bins are split when using manual binning
- `adaptive`: determines if the R23 bins have equal or different number of spectra in them in binning method
- `dustatten`: determines if dust attenuation corrections are applied
- `mask`: determines if the night sky mask is used in Stackingboth_MasterGrid.py

Difference between analyses

Different grid methods were utilize throughout the process of developing this study. The dataset option determines
which grid method is used. For the run_grid_R23_O32_analysis(), the following options are available. 
- `Grid`: set two dimensional grid with equal lengthen sides 
- `O32_Grid`: one dimensional grid that bins in O32 
- `R23_Grid`: one dimensional grid that bins in R23
- `n_Bins`: two dimensional grid with a set number of spectra in each bin
            bins in R23 then in O32

Sample parameters for general.run_grid_r23_o32_analysis()
```python
dataset = 'n_Bins'
y_correction = ''
n_split = 3
adaptive = True 
dustatten = True
mask = True
```

Calling the run function:

```python 
    general.run_grid_r23_o32_analysis(dataset, y_correction, n_split, 
                                                               adaptive, dustatten, mask)
```
    
Steps taking throughout run function: 

1. Gets the valid data for the study using get_det3() 
```python 
    general.get_det3(fitspath, fitspath_ini)
```
    
2. Calls correct binning function for dataset 
    
    Example: 
```python
    if dataset == 'n_Bins':
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_pages, bin_outfile, n_split, individual_ID,
                                            R23, O32, SNR3, data3, galinbin)
```

3. Calls stacking function to stack individual spectra 
    
    A mask can be applied to correct for night sky lines. Produces a table with binned data properties. 

```python 
    stackboth_mastergrid.run_stacking_master_mask(fitspath, fitspath_ini, dataset, stack_name, bin_outfile)
```

4. Calls fitting function to fit the emisison lines in the combined spectra 
with a gaussian profile and determine gaussian properties 
   
   Produces a table with emission line fitting parameters for each line of each binned spectra. 
```python 
    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave, lineflag, dispersion, y_correction,
                                      tab=binning_avg_asc)
```

5. Creates a validation table that confirms detections of OIIIλ[4363]

   Code imported from MSC. Adds a row to the table of emission line measurements to indicate if there is a detection of OIIIλ[4363].
    A one represents a detection, while a zero represents a nan-detection. 
```python 
    # Verification Table
    valid_table.make_validation_table(fitspath)
```
    
6. Calls function to calculate the R value, temperature, and metallicity of the detected lines
   
   Calculates the _R_ flux ratio, electron temperature, and metallicities of each bin and saves in table. 
```python 
    r_temp_calcul.run_function(fitspath, dataset, verification_table_revised, dustatt=False)
```

7. Applies dust attenuation corrections if specified in run function 
   
```python 
    if dustatten:
        balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits', use_revised=False)
        attenuation.EBV_table_update(fitspath, use_revised= False)
        r_temp_calcul.run_function(fitspath, dataset, verification_table_revised, dustatt=True)
```

8. Applies Error Propagation from MSC
   
```python 
    error_prop.fluxes_derived_prop(fitspath, binned_data=True, revised=True)
```

9. Creates calibration plots that compare detections to other studies 
   
```python 
    calibration_plots.lac_gpc_plots(fitspath, fitspath_ini, dataset, revised=True, individual=False)
```

### Running Voronoi Analysis 
The analysis of the binning methods is run by executing the `run_grid_R23_O32_analysis()`
function in [analysis/deep2_r23_o32/archives/run_functions/voronoi_general.py](deep2_r23_o32/archives/run_functions/voronoi_general.py). 
```python
    from Zcalbase_gal.analysis.deep2_r23_o32.archives.run_functions import voronoi_general
```

The `run` function requires the following inputs:

- `dataset`: keyword used to define binning method options: 'Voronoi10', 'Voronoi14', 'Voronoi20'
- `y_correction`: determines if the smoothed (movingaverage_box1D) version of y is used in zoom_and_gauss_general.py
- `dustatten`: determines if dust attenuation corrections are applied
- `mask`: determines if the night sky mask is used in Stackingboth_MasterGrid.py

Difference between analyses
Different grid methods were utilize throughout the process of developing this study including the Voronoi 
Tessellation code. The dataset option for voronoi run function determines the target single to noise of each bin,
which varies the number of spectra in each bin. 

Sample parameters for `general.run_grid_r23_o32_analysis()`
```python
dataset = 'Voronoi14'
y_correction = ''
dustatten = True
mask = True
```

Calling the run function:
```python
    general.run_voronoi_r23_o32_analysis(dataset, y_correction, n_split, adaptive, dustatten, mask)
```
This run function goes through the same process as the grid method above. 
