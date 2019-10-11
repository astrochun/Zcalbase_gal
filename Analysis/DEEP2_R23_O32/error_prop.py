

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
def error_prop_chuncodes(fitspath, project, dataset):
    if project == 'R23O32':
        TM_tab = asc.read(fitspath +dataset +'_temperatures_metalicity.tbl')
        combine_flux_tab = asc.read(fitspath +dataset+'_combined_flux_table.tbl')
        
        #R23O32 table calls
        OII_flux      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
        OII_RMS       = combine_flux_tab['OII_3727_RMS'].data
        NeIII_flux    = combine_flux_tab['NeIII_Flux_Gaussian'].data
        NeIII_RMS     = combine_flux_tab['NeIII_RMS'].data
        HeI_flux      = combine_flux_tab['HeI_Flux_Gaussian'].data
        HeI_RMS       = combine_flux_tab['HeI_RMS'].data
        Hdelta_flux   = combine_flux_tab['HDELTA_Flux_Gaussian'].data
        Hdelta_RMS    = combine_flux_tab['HDELTA_RMS'].data
        Hgamma_flux   = combine_flux_tab['Hgamma_Flux_Gaussian'].data
        Hgamma_RMS    = combine_flux_tab['Hgamma_RMS'].data
        OIII4363_flux = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
        OIII4363_RMS  = combine_flux_tab['OIII_4363_RMS'].data
        OIII4959_flux = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
        OIII4959_RMS  = combine_flux_tab['OIII_4958_RMS'].data
        OIII5007_flux = combine_flux_tab['OIII_5007_Flux_Gaussian'].data
        OIII5007_RMS  = combine_flux_tab['OIII_5007_RMS'].data
        
        flux_data = [OII_flux,NeIII_flux, HeI_flux, Hdelta_flux, Hgamma_flux,OIII4363_flux,OIII4959_flux, OIII5007_flux ]
        RMS_data  = [OII_RMS,NeIII_RMS, HeI_RMS, Hdelta_RMS, Hgamma_RMS,OIII4363_RMS,OIII4959_RMS, OIII5007_RMS]
        
    if project == 'ML':
        TM_tab = asc.read(fitspath + dataset + '_derived_properties_metallicity.tbl')
        combine_flux_tab = asc.read(fitspath + dataset + '_emission_lines.tbl')
        
        print('ML table calls')
        #ML table calls
        OII_flux      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
        OII_RMS       = combine_flux_tab['OII_3727_RMS'].data
        Hdelta_flux   = combine_flux_tab['HDELTA_Flux_Gaussian'].data
        Hdelta_RMS    = combine_flux_tab['HDELTA_RMS'].data
        Hgamma_flux   = combine_flux_tab['HGAMMA_Flux_Gaussian'].data
        Hgamma_RMS    = combine_flux_tab['HGAMMA_RMS'].data
        OIII4363_flux = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
        OIII4363_RMS  = combine_flux_tab['OIII_4363_RMS'].data
        OIII4959_flux = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
        OIII4959_RMS  = combine_flux_tab['OIII_4958_RMS'].data
        OIII5007_flux = combine_flux_tab['OIII_5007_Flux_Gaussian'].data
        OIII5007_RMS  = combine_flux_tab['OIII_5007_RMS'].data
        
        flux_data = [OII_flux, Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux ]
        RMS_data  = [OII_RMS, Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]
        
        
    if project == '': print('Please specify project')
    

    Temp          = TM_tab['Temperature'].data
    metallicity   = TM_tab['com_O_log'].data

    x_arr = np.random.random(1000)
    #p_page= fitspath+dataset+'/'+str(np.int(working_wave))+'error_histogram.pdf'
    
    #pdf_pages_err = PdfPages(p_page)
    
    for aa in range(len(flux_data)):
        flux_gpdf = random_pdf(flux_data[aa],RMS_data[aa], seed_i =1, n_iter=1000, silent = False)

    print('err_function:', flux_gpdf, flux_gpdf.shape)

    
    for bb in range(len(flux_data)): 
        err, xpeak= compute_onesig_pdf(flux_gpdf,flux_data[bb], usepeak=False, silent=True, verbose = True)
    print('err',err, len(err),'xpeak', xpeak,len(err))

    ###Next Step is to make the curves and graph the functions

    
    #fig_err, ax = plt.subplots()
    #plt.hist(fluxg_pdf, 50)

    #fig_err.set_size_inches(8,8)
    #fig_err.savefig(pdf_pages_err, format='pdf')
    #pdf_pages_err.close()
    
