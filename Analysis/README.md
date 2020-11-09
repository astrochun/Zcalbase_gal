# Executing `Zcalbase_gal` on DEEP2 data

1. [Overview](#overview)
2. [Execution](#execution)
    1. [Requirements](#requirements)
    2. [Running Grid Analysis](#running-grid-analysis)
    3. [Running Voronoi Analysis](#running-voronoi-analysis)

# Overview

Python 3.x codes analysis of a subset of the [DEEP2 Redshift Survey](http://deep.ps.uci.edu).
Code relies on [Metallicity Stack Commons (MSC)](https://github.com/astrochun/Metallicity_Stack_Commons)
created by Chun Ly, Reagen Leimbach, and Caroline McCormick. 

Previous studies of galaxy evolution have been proven incorrect due to the fact
that the study did not take into account the difference in the gas properties
between high- and low-z galaxies (redshifted). This analysis presents a
calibration between the metallicity and strong-line diagnostics, R23 and O32,
to use to correct for the change in gas properties as redshift increases. To
create this calibration, the OIIIλ[4363] is measured to calculate the electron
temperature of the galaxy, which is used to calculate the metallicity.

OIIIλ[4363] is a forbidden transition, which produces a very weak emission line
in the spectra. In order to detect OIIIλ[4363], individual spectra are binned
and stacked to get a composite spectrum measurement.

This analysis is also compared against recently published analyses from:
1. [Jiang, T., Malhotra, S., Rhoads, J. E., et al. 2019, ApJ, 872, 145](https://arxiv.org/abs/1811.05796)  
2. [Bian, F., Kewley, L. J., & Dopita, M. A. 2018, ApJ, 859, 175](https://iopscience.iop.org/article/10.3847/1538-4357/aabd74/meta)

Links to material presenting this work: 
1. [UA Space Grant Link](https://arizona.figshare.com/articles/Stacking_of_Galaxy_Spectra/12360626) 
to figure created in the stacking code 


2. [Link (TBD)] to Honors Thesis

# Execution 

### Requirements 
Your will need the following to have a working copy of this software.

- Python (>=3.7). We used 3.7.7
- [Metallicity Stack Commons](https://github.com/astrochun/Metallicity_Stack_Commons)
- numpy 
- matplotlib
- astropy
- scipy
- [chun_codes](https://github.com/astrochun/chun_codes)

### Running Grid Analysis 
The analysis of the binning methods is run by executing the
`run_grid_R23_O32_analysis` function in
[`analysis.deep2_r23_o32.general`](deep2_r23_o32/general.py).

``` python
from Zcalbase_gal.analysis.deep2_r23_o32 import general
```

The `run_grid_R23_O32_analysis` function requires the following variables:
- `dataset`: variable used to define binning method options:
  'Grid', 'O32_Grid', 'R23_Grid', or 'n_Bins'. See [next section](#different-grid-methods)
- `y_correction`: Bool. Determines if the smoothed (`movingaverage_box1D`)
  version of spectra is used in `zoom_and_gauss_general` module
- `n_split`: determined how many times the log(R23) bins are split in log(O32)
  when using manual binning
- `adaptive`: Bool. Set for log(R23) bins to have equal number of spectra in
  them in binning method. Default: `False`
- `dustatten`: Bool. Set to apply dust attenuation corrections. Default: `False`
- `mask`: Indicate if night sky masking is implemented in
  `Stackingboth_MasterGrid`. Default: `None`/`False`

#### Different grid methods

Different grid methods were utilize throughout the process of developing this
study. The `dataset` option determines which grid method is used. For the
`run_grid_R23_O32_analysis`, the following options are available:
 - `Grid`: Two-dimensional grid with equal bin size
 - `O32_Grid`: One-dimensional grid that bins in log(O32)
 - `R23_Grid`: One-dimensional grid that bins in log(R23)
 - `n_Bins`: Two-dimensional grid with a set number of spectra in each
             log(R23) bin followed by sub-binning in log(O32)

Sample parameters for `general.run_grid_r23_o32_analysis`
``` python
dataset = 'n_Bins'
y_correction = False
n_split = 3
adaptive = True 
dustatten = True
mask = True
```

Calling the `run_grid_R23_O32_analysis` function:

``` python
general.run_grid_r23_o32_analysis(dataset, y_correction, n_split,
                                  adaptive, dustatten, mask)
```
    
Steps taking throughout the `run_grid_R23_O32_analysis` function:

1. Gets the valid data for the study using `get_det3` function:

   ``` python
   general.get_det3(fitspath, fitspath_ini)
   ```
    
2. Calls the correct binning function for dataset

   For example:
   ``` python
   if dataset == 'n_Bins':
      n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_pages, bin_outfile,
                                          n_split, individual_ID, R23, O32,
                                          SNR3, data3, galinbin)
   ```

3. Calls the stacking function to stack individual spectra
    
   This produces a table with binned data properties.
   A mask can be applied to correct for night sky lines. 

   ``` python
   # With masking
   stackboth_mastergrid.run_stacking_master_mask(fitspath, fitspath_ini, dataset,
                                                 name, grid_data_file)

   # Without masking
   stackboth_mastergrid.run_stacking_master(fitspath, name, grid_data_file):
   ```

4. Calls fitting function to fit the emission lines in the combined spectra
with a gaussian profile and determine Gaussian fitting parameters
 
   This generates a table with emission-line fitting results for each line of
   each binned spectra.

   ``` python
   zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave, lineflag,
                                     dispersion, y_correction, tab=binning_avg_asc)
   ```

5. Creates a validation table that confirms detections of OIIIλ[4363]

   Code is imported from MSC. It creates a column to the table of emission line
   measurements to indicate if there is a detection of OIIIλ[4363].
   Detection: 1; Non-detection: 0
   ``` python
   from Metallicity_Stack_commons import valid_table
   valid_table.make_validation_table(fitspath)
   ```
    
6. Calls function to calculate the _R_ value, temperature, and metallicity of
   the detected lines
   
   Calculates the _R_ flux ratio, electron temperature, and metallicities of
   each bin and saves in table.
   ``` python
   r_temp_calcul.run_function(fitspath, dataset, verification_table_revised,
                              dustatt=False)
   ```

7. Applies dust attenuation corrections if specified in run function 
   
   ``` python
   from Metallicity_Stack_commons.plotting import balmer
   if dustatten:
       balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits', use_revised=False)
       attenuation.EBV_table_update(fitspath, use_revised= False)
       r_temp_calcul.run_function(fitspath, dataset, verification_table_revised,
                                  dustatt=True)
   ```

8. Applies error propagation
   
   ``` python
   from Metallicity_Stack_commons.analysis import error_prop
   error_prop.fluxes_derived_prop(fitspath, binned_data=True, revised=True)
   ```

9. Creates calibration plots that compare detections to other studies 
   
   ``` python
   calibration_plots.lac_gpc_plots(fitspath, fitspath_ini, dataset,
                                   revised=True, individual=False)
   ```

### Running Voronoi Analysis 
The analysis of the binning methods is run by executing the `run_grid_R23_O32_analysis`
function in [`analysis.deep2_r23_o32.archives.run_functions.voronoi_general`](deep2_r23_o32/archives/run_functions/voronoi_general.py).
``` python
from Zcalbase_gal.analysis.deep2_r23_o32.archives.run_functions import voronoi_general
```

The `run_grid_R23_O32_analysis` function requires the following inputs:

- `dataset`: Variable used to define binning method options: 'Voronoi10', 'Voronoi14', 'Voronoi20'
- `y_correction`: determines if the smoothed (movingaverage_box1D) version of y is used in zoom_and_gauss_general.py
- `dustatten`: determines if dust attenuation corrections are applied
- `mask`: determines if the night sky mask is used in Stackingboth_MasterGrid.py

#### Difference between analyses

Different grid methods were utilize throughout the process of developing this
study including the Voronoi Tessellation code. The dataset option for Voronoi
run function determines the target single-to-noise of each bin, which varies
the number of spectra in each bin.

Sample parameters for `general.run_grid_r23_o32_analysis`
``` python
dataset = 'Voronoi14'
y_correction = False
dustatten = True
mask = True
```

Calling the `run_grid_R23_O32_analysis` function:
``` python
general.run_voronoi_r23_o32_analysis(dataset, y_correction, n_split, adaptive,
                                     dustatten, mask)
```
This run function goes through the same process as the grid method above. 
