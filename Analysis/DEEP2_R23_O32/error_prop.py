

import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf


#line_name = ['OII_3727','NeIII','HeI', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']

#Calling Error Propogation
#error_prop_chuncodes(fitspath, dataset, working_wave, flux_g_array, RMS_array)
# TM_tab is the table with your metallicities, temperatures, and line ratios
# flux_tab (for Reagen) is produced by the zoom_and_gauss_general code and has the RMS values for each bin
# flux_tab (for Caroline) is produced by emission_line_fit code and has the RMS values for each bin
def error_prop_chuncodes(fitspath, flux_file):
    #TM_file = fitspath + dataset + '_temperatures_metalicity.tbl'
    #flux_file = fitspath + dataset + '_combined_flux_table.tbl'
    #TM_file = fitspath + dataset + '_derived_properties_metallicity.tbl'
    #flux_file = fitspath + dataset + '_emission_lines.tbl'
    
    combine_flux_tab = asc.read(flux_file)
    
    OII_flux      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
    OII_RMS       = combine_flux_tab['OII_3727_RMS'].data
    Hdelta_flux   = combine_flux_tab['HDELTA_Flux_Gaussian'].data
    Hdelta_RMS    = combine_flux_tab['HDELTA_RMS'].data
    OIII4363_flux = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
    OIII4363_RMS  = combine_flux_tab['OIII_4363_RMS'].data
    OIII4959_flux = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
    OIII4959_RMS  = combine_flux_tab['OIII_4958_RMS'].data
    OIII5007_flux = combine_flux_tab['OIII_5007_Flux_Gaussian'].data
    OIII5007_RMS  = combine_flux_tab['OIII_5007_RMS'].data
    Hgamma_flux   = combine_flux_tab['HGAMMA_Flux_Gaussian'].data
    Hgamma_RMS    = combine_flux_tab['HGAMMA_RMS'].data
        
    
    line_names = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'OIII_4958', 'OIII_5007']
    flux_data = [OII_flux, Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux ]
    RMS_data  = [OII_RMS, Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]

    x_arr = np.random.random(1000)
    #p_page= fitspath+dataset+'/'+str(np.int(working_wave))+'error_histogram.pdf'
    
    #pdf_pages_err = PdfPages(p_page)
    
    for aa in range(len(flux_data)):
        flux_gpdf = random_pdf(flux_data[aa],RMS_data[aa], seed_i =aa, n_iter=1000, silent = False)
        err, xpeak= compute_onesig_pdf(flux_gpdf,flux_data[aa], usepeak=False, silent=True, verbose = True)
        
        ##Won't work for Reagen's Hgamma_Flux_Gaussian column bc lower case gamma
        col_name_idx = combine_flux_tab.index_column(line_names[aa] + '_Flux_Gaussian')
        err_t = err.transpose()
        c1 = Column(err_t[0], name=line_names[aa]+'_Low_Error')
        c2 = Column(err_t[1], name=line_names[aa]+'_High_Error')
        combine_flux_tab.add_columns([c1, c2], indexes=[col_name_idx,col_name_idx])
   
        print('err_function:', flux_gpdf, flux_gpdf.shape)
        print('err',err, len(err),'xpeak', xpeak,len(err))
    asc.write(combine_flux_tab, flux_file, format= 'fixed_width_two_line')
        

    ###Next Step is to make the curves and graph the functions
    
def plotting_errors(fitspath, TM_file, flux_file):
       
    TM_tab = asc.read(TM_file)
    combine_flux_tab = asc.read(flux_file)

    Temp          = TM_tab['Temperature'].data
    metallicity   = TM_tab['com_O_log'].data

    OII_flux      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
    OIII4363_flux = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
    OIII4959_flux = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
    OIII5007_flux = combine_flux_tab['OIII_5007_Flux_Gaussian'].data

    OII_Low_Error      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
    OIII4363_Low_Error = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
    OIII4959_Low_Error = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
    OIII5007_Low_Error = combine_flux_tab['OIII_5007_Flux_Gaussian'].data

    OII_High_Error      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
    OIII4363_High_Error = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
    OIII4959_High_Error = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
    OIII5007_High_Error = combine_flux_tab['OIII_5007_Flux_Gaussian'].data
    
    
    


    #fig_err, ax = plt.subplots()
    #plt.hist(fluxg_pdf, 50)

    #fig_err.set_size_inches(8,8)
    #fig_err.savefig(pdf_pages_err, format='pdf')
    #pdf_pages_err.close()
    
