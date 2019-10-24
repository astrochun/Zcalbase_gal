

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

#Imports R_temp_calcul code for temp/metallicity calculation
from Zcalbase_gal.Analysis.DEEP2_R23_O32 import R_temp_calcul


#line_name = ['OII_3727','NeIII','HeI', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']

#Calling Error Propogation
#error_prop_chuncodes(fitspath, dataset, working_wave, flux_g_array, RMS_array)
# TM_tab is the table with your metallicities, temperatures, and line ratios
# flux_tab (for Reagen) is produced by the zoom_and_gauss_general code and has the RMS values for each bin
# flux_tab (for Caroline) is produced by emission_line_fit code and has the RMS values for each bin
def error_prop_chuncodes(fitspath, flux_file, TM_file):
    #TM_file = fitspath + dataset + '_temperatures_metalicity.tbl'
    #flux_file = fitspath + dataset + '_combined_flux_table.tbl'
    #TM_file = fitspath + dataset + '_derived_properties_metallicity.tbl'
    #flux_file = fitspath + dataset + '_emission_lines.tbl'
    
    combine_flux_tab = asc.read(flux_file)
    TM_tab = asc.read(TM_file)
    
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
    Hbeta_flux    = combine_flux_tab['HBETA_Flux_Gaussian'].data
    Hbeta_RMS     = combine_flux_tab['HBETA_RMS'].data

    Temp          = TM_tab['Temperature'].data
    com_O_log     = TM_tab['com_O_log'].data
    O_s_ion       = TM_tab['O_s_ion'].data
    O_d_ion       = TM_tab['O_d_ion'].data
    log_O_s       = TM_tab['log_O_s'].data
    log_O_d       = TM_tab['log_O_d'].data

    
    line_names = ['OII_3727', 'HBETA', 'HDELTA', 'HGAMMA', 'OIII_4363', 'OIII_4958', 'OIII_5007']
    flux_data = [OII_flux,Hbeta_flux ,Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux ]
    RMS_data  = [OII_RMS, Hbeta_RMS,Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]

    x_arr = np.random.random(1000)
    #p_page= fitspath+dataset+'/'+str(np.int(working_wave))+'error_histogram.pdf'

    #Initialize Dictionary for flux_gpdf
    pdf_dict = {}
    
    for aa in range(len(flux_data)):
        flux_gpdf = random_pdf(flux_data[aa],RMS_data[aa], seed_i =aa, n_iter=1000, silent = False)
        err, xpeak= compute_onesig_pdf(flux_gpdf,flux_data[aa], usepeak=False, silent=True, verbose = True)

        #Fill In Dictionary
        pdf_dict[line_names[aa]] = flux_gpdf
        
        '''col_name_idx = combine_flux_tab.index_column(line_names[aa] + '_Flux_Gaussian')
        err_t = err.transpose()
        c1 = Column(err_t[0], name=line_names[aa]+'_Low_Error')
        c2 = Column(err_t[1], name=line_names[aa]+'_High_Error')
        combine_flux_tab.add_columns([c1, c2], indexes=[col_name_idx,col_name_idx])'''
   
        print('err_function:', flux_gpdf, flux_gpdf.shape)
        print('err',err, len(err),'xpeak', xpeak,len(err))
    asc.write(combine_flux_tab, flux_file, format= 'fixed_width_two_line')
    print pdf_dict

    ####################R_temp_calcul calls############################
    EBV = np.zeros(pdf_dict['OIII_4363'].shape)
    k_4363 = np.zeros(pdf_dict['OIII_4363'].shape)
    k_5007 = np.zeros(pdf_dict['OIII_4363'].shape)

    R_pdf = R_temp_calcul.R_calculation(pdf_dict['OIII_4363'], pdf_dict['OIII_5007'], pdf_dict['OIII_4958'], EBV, k_4363, k_5007)

    Te_pdf = R_temp_calcul.temp_calculation(R_pdf)

    two_beta_pdf = pdf_dict['OII_3727']/pdf_dict['HBETA']
    three_beta_pdf = (pdf_dict['OIII_5007'] + pdf_dict['OIII_4958'])/pdf_dict['HBETA']

    O_s_ion_pdf , O_d_ion_pdf, com_O_log_pdf, O_s_ion_log_pdf, O_d_ion_log_pdf = R_temp_calcul.metalicity_calculation(Te_pdf, two_beta_pdf, three_beta_pdf)

    print O_s_ion.shape

    metal_error = {}
    metal_xpeak = {}
    
    ###########compute_onesig_pdf for all the metallicity outputs#############
    # Pass in the pdf metallicities, all the stacked measurements
    metal_str = ['O_s_ion_pdf' , 'O_d_ion_pdf', 'com_O_log_pdf', 'O_s_ion_log_pdf', 'O_d_ion_log_pdf']
    metallicity_names = [O_s_ion_pdf , O_d_ion_pdf, com_O_log_pdf, O_s_ion_log_pdf, O_d_ion_log_pdf]
    
    combined_metallicity = [O_s_ion , O_d_ion, com_O_log, log_O_s, log_O_d]
    for ii in range(len(metallicity_names)):
        err_metal, xpeak_metal = compute_onesig_pdf(metallicity_names[ii], combined_metallicity[ii], usepeak=False, silent=True, verbose = True)
        print err_metal, len(err_metal), xpeak_metal

        metal_error[metal_str[ii]] = err_metal
        metal_xpeak[metal_str[ii]] = xpeak_metal

    print metal_error
    print metal_xpeak


    np.savez(fitspath+'metal_errors.npz', **metal_error)
    np.savez(fitspath+'metal_xpeaks.npz', **metal_xpeak)
    

    
def plotting_errors(fitspath, flux_file):
    pdf_pages = PdfPages(fitspath+'emission_line_error_graphs.pdf')
    #TM_tab = asc.read(TM_file)
    combine_flux_tab = asc.read(flux_file)

    
    
    ID = combine_flux_tab['ID'].data
    OII_flux      = combine_flux_tab['OII_3727_Flux_Gaussian'].data
    OIII4363_flux = combine_flux_tab['OIII_4363_Flux_Gaussian'].data
    OIII4959_flux = combine_flux_tab['OIII_4958_Flux_Gaussian'].data
    OIII5007_flux = combine_flux_tab['OIII_5007_Flux_Gaussian'].data

    OII_Low_Error      = np.array([combine_flux_tab['OII_3727_Low_Error'].data])
    OIII4363_Low_Error = np.array([combine_flux_tab['OIII_4363_Low_Error'].data])
    OIII4959_Low_Error = np.array([combine_flux_tab['OIII_4958_Low_Error'].data])
    OIII5007_Low_Error = np.array([combine_flux_tab['OIII_5007_Low_Error'].data])

    OII_High_Error      = np.array([combine_flux_tab['OII_3727_High_Error'].data])
    OIII4363_High_Error = np.array([combine_flux_tab['OIII_4363_High_Error'].data])
    OIII4959_High_Error = np.array([combine_flux_tab['OIII_4958_High_Error'].data])
    OIII5007_High_Error = np.array([combine_flux_tab['OIII_5007_High_Error'].data])
    
    error_OII = np.array([OII_Low_Error, OII_High_Error], dtype=np.float)
    error_OIII4363 = np.array([OIII4363_Low_Error, OIII4363_High_Error], dtype=np.float)
    error_OIII4959 = np.array([OIII4959_Low_Error, OIII4959_High_Error], dtype=np.float)
    error_OIII5007 = np.array([OIII5007_Low_Error, OIII5007_High_Error], dtype=np.float)

    
    line_names = ['OII_3727', 'OIII_4363', 'OIII_4958', 'OIII_5007']

    
    
    for ii in range(len(ID)):
        ##plot Gaussian
        fig,ax = plt.subplots()
        ax.plot(3727, OII_flux[ii], 'r', marker= 'o', markersize = 2)
        ax.plot(4363, OIII4363_flux[ii], 'r', marker= 'o', markersize = 2)
        ax.plot(4959, OIII4959_flux[ii], 'r', marker= 'o', markersize = 2)
        ax.plot(5007, OIII5007_flux[ii], 'r', marker= 'o', markersize = 2)
        
        ax.errorbar(3727, OII_flux[ii], fmt= 'None', yerr = error_OII, ecolor = 'k')
        ax.errorbar(4363, OIII4363_flux[ii], fmt= 'None',yerr = error_OIII4363, ecolor = 'k')
        ax.errorbar(4959, OIII4959_flux[ii],fmt= 'None', yerr = error_OIII4959, ecolor = 'k')
        ax.errorbar(5007, OIII5007_flux[ii], fmt= 'None',yerr = error_OIII5007, ecolor = 'k')        
        

        fig.savefig(pdf_pages, format='pdf')
        
    fig.set_size_inches(8,8)

    
    pdf_pages.close()
    
    
    


    #fig_err, ax = plt.subplots()
    #plt.hist(fluxg_pdf, 50)

    #fig_err.set_size_inches(8,8)
    #fig_err.savefig(pdf_pages_err, format='pdf')
    #pdf_pages_err.close()
    
