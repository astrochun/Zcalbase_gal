####Runs ranom_pudf and compute_onesig_pdf from chun_codes to perform error propagation######
####Error propagation for fluxes, temperature, and metallicity#####
'''

################Functions###################
#error_prop_chuncodes(path, flux_file, TM_file)
Input variables: 
path           -> location of where the outputted pdf_file will be saved
flux_file      -> ascii file with all the line fluxes and RMS values created by 
                  zoom_and_gauss_general (Reagen) and emission_line_fit (Caroline)
TM_file        -> ascii file with the temperature and metallicity of the bins 
                  created by R_temp_calcul.py

Outputs: Four dictionaries of the temperatures and metallicities 

#plotting_linear_errors(fitspath, flux_file)
Function not being called currently
Output: pdf file with plots of calculated fluxes with lower and high bounds marked 
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
from pylab import subplots_adjust
from scipy.optimize import curve_fit
import scipy.integrate as integ

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf

#Imports R_temp_calcul code for temp/metallicity calculation
#from Zcalbase_gal.Analysis.DEEP2_R23_O32 import R_temp_calcul
from Metallicity_Stack_Commons.temp_metallicity_calc import R_calculation, temp_calculation, metallicity_calculation

from Metallicity_Stack_Commons.column_names import temp_metal_names0, bin_names0


#line_name = ['OII_3727','NeIII','HeI', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']

def error_prop_chuncodes(path, flux_file, TM_file, verification_tab):    
    #TM_file = fitspath + dataset + '_temperatures_metalicity.tbl'
    #flux_file = fitspath + dataset + '_combined_flux_table.tbl'
    #TM_file = fitspath + dataset + '_derived_properties_metallicity.tbl'
    #flux_file = fitspath + dataset + '_emission_lines.tbl'
    
    combine_flux_tab = asc.read(flux_file)
    TM_tab = asc.read(TM_file)
    verify = asc.read(verification_tab)

    detect = verify[bin_names0[2]]
    ID = verify[bin_names0[0]]
    detection = np.where((detect==1))[0]
    ID_detect = ID[detection]
    print (ID_detect)

    
    combine_flux = combine_flux_tab[detection]
    TM = TM_tab[detection]

    OII_flux      = combine_flux['OII_3727_Flux_Gaussian'].data
    OII_RMS       = combine_flux['OII_3727_RMS'].data
    Hdelta_flux   = combine_flux['HDELTA_Flux_Gaussian'].data
    Hdelta_RMS    = combine_flux['HDELTA_RMS'].data
    OIII4363_flux = combine_flux['OIII_4363_Flux_Gaussian'].data
    OIII4363_RMS  = combine_flux['OIII_4363_RMS'].data
    OIII4959_flux = combine_flux['OIII_4958_Flux_Gaussian'].data
    OIII4959_RMS  = combine_flux['OIII_4958_RMS'].data
    OIII5007_flux = combine_flux['OIII_5007_Flux_Gaussian'].data
    OIII5007_RMS  = combine_flux['OIII_5007_RMS'].data
    Hgamma_flux   = combine_flux['HGAMMA_Flux_Gaussian'].data
    Hgamma_RMS    = combine_flux['HGAMMA_RMS'].data
    Hbeta_flux    = combine_flux['HBETA_Flux_Gaussian'].data
    Hbeta_RMS     = combine_flux['HBETA_RMS'].data


    Temp          = TM[temp_metal_names0[0]].data
    com_O_log     = TM[temp_metal_names0[1]].data
    O_s_ion       = TM[temp_metal_names0[4]].data
    O_d_ion       = TM[temp_metal_names0[5]].data
    log_O_s       = TM[temp_metal_names0[2]].data
    log_O_d       = TM[temp_metal_names0[3]].data

      
    #################Starting Error Calculations#################################################
    line_names = ['OII_3727', 'HBETA', 'HDELTA', 'HGAMMA', 'OIII_4363', 'OIII_4958', 'OIII_5007']
    flux_data = [OII_flux,Hbeta_flux ,Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux ]
    RMS_data  = [OII_RMS, Hbeta_RMS,Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]

    #p_page= fitspath+dataset+'/'+str(np.int(working_wave))+'error_histogram.pdf'

    #Initialize Dictionary for flux_gpdf
    flux_propdist_dict = {}
    flux_xpeak = {}
    flux_lowhigh = {}
    
    for aa in range(len(flux_data)):
        flux_gpdf = random_pdf(flux_data[aa],RMS_data[aa], seed_i =aa, n_iter=1000, silent = False)
        err, xpeak= compute_onesig_pdf(flux_gpdf,flux_data[aa], usepeak=True, silent=True, verbose = True)

        #Fill In Dictionary
        flux_propdist_dict[line_names[aa]+'_pdf'] = flux_gpdf
        flux_xpeak[line_names[aa]+'_xpeak'] = xpeak
        flux_lowhigh[line_names[aa]+'_lowhigh_error'] = err

        #Edit Ascii Table
        combine_flux_new_name = flux_file.replace('.tbl','.revised.tbl')
        
        combine_flux_tab[line_names[aa]+'_Flux_Gaussian'][detection] = xpeak 
        
        #col_name_idx = combine_flux_tab.index_column(line_names[aa] + '_Flux_Gaussian')
        #err_t = err.transpose()
        #c1 = Column(err_t[0], name=line_names[aa]+'_Low_Error')
        #c2 = Column(err_t[1], name=line_names[aa]+'_High_Error')
        #combine_flux_tab.add_columns([c1, c2], indexes=[col_name_idx,col_name_idx])
   
        #print('err_function:', flux_gpdf, flux_gpdf.shape)
        #print('err',err, len(err),'xpeak', xpeak,len(err))
    asc.write(combine_flux_tab, combine_flux_new_name, format= 'fixed_width_two_line')
    #print pdf_dict

    ####################R_temp_calcul calls with probability distribution of fluxes############################
    EBV = np.zeros(flux_propdist_dict['OIII_4363_pdf'].shape)
    k_4363 = np.zeros(flux_propdist_dict['OIII_4363_pdf'].shape)
    k_5007 = np.zeros(flux_propdist_dict['OIII_4363_pdf'].shape)

    Te_propdist_dict = {}
    Te_xpeaks = {}
    Te_lowhigh = {}
    #changed parameters for R_calculation --> got rid of 4959

    R_propdist = R_calculation(flux_propdist_dict['OIII_4363_pdf'], flux_propdist_dict['OIII_5007_pdf'], EBV)

    Te_propdist = temp_calculation(R_propdist)
    err_te, xpeak_te = compute_onesig_pdf(Te_propdist, Temp, usepeak=True, silent=True, verbose = True)

    Te_propdist_dict['Te_pdf'] = Te_propdist
    Te_xpeaks['Te_xpeak']= xpeak_te
    Te_lowhigh['Te_lowhigh_error'] = err_te

    two_beta = flux_propdist_dict['OII_3727_pdf']/flux_propdist_dict['HBETA_pdf']
    #changed three_beta to 5007(1 + 1/3.1) instead of 5007 + 4959
    three_beta = (flux_propdist_dict['OIII_5007_pdf'] * (1 + 1/3.1))/flux_propdist_dict['HBETA_pdf']

    #Edit TM_File to add temperature errors
    TM_new_name = TM_file.replace('.tbl','.revised.tbl')
    TM_tab[temp_metal_names0[0]][detection] = xpeak_te
    
    #err_tT = err_te.transpose()
    #c1 = Column(err_tT[0], name='Temperature_Low_Error')
    #c2 = Column(err_tT[1], name='Temperature_High_Error')
    #TM_tab.add_columns([c1, c2], indexes=[18,18])
    
    com_O_log_propdist, metal_dict = metallicity_calculation(Te_propdist, two_beta, three_beta)

    O_s_ion_log_propdist = metal_dict['log(O+/H)']  
    O_d_ion_log_propdist = metal_dict['log(O++/H)']
    O_s_ion_propdist = metal_dict['O+/H']
    O_d_ion_propdist = metal_dict['O++/H']

    com_O_log_propdist, metal_dict = R_temp_calcul.metallicity_calculation(Te_propdist, two_beta, three_beta)
    metallicity_propdist_dict = {'O_s_ion_pdf': metal_dict['O_s_ion'], 'O_d_ion_pdf': metal_dict['O_d_ion'] , 'com_O_log_pdf': com_O_log_propdist, 'O_s_ion_log_pdf': metal_dict['O_s_ion_log'], 'O_d_ion_log_pdf':  metal_dict['O_d_ion_log']}
    
    
    metal_error = {}
    metal_xpeak = {}
    
    
    ###########compute_onesig_pdf for all the metallicity outputs#############
    # Pass in the pdf metallicities, all the stacked measurements
    metal_str = ['O_s_ion_pdf' , 'O_d_ion_pdf', 'com_O_log_pdf', 'O_s_ion_log_pdf', 'O_d_ion_log_pdf']
    metal_err = ['O_s_ion_lowhigh_error' , 'O_d_ion_lowhigh_error', 'com_O_log_lowhigh_error', 'O_s_ion_log_lowhigh_error', 'O_d_ion_log_lowhigh_error']
    xpeak_key = ['O_s_ion_xpeak','O_d_ion_xpeak', 'com_O_log_xpeak', 'O_s_ion_log_xpeak', 'O_d_ion_log_xpeak']

    #metallicity_names =[O_s_ion_propdist , O_d_ion_propdist, com_O_log_propdist, O_s_ion_log_propdist, O_d_ion_log_propdist]
    metallicity_names= [metal_dict['O_s_ion'], metal_dict['O_d_ion'] ,  com_O_log_propdist, metal_dict['O_s_ion_log'], metal_dict['O_d_ion_log']]
    
    combined_metallicity = [O_s_ion , O_d_ion, com_O_log, log_O_s, log_O_d]
    metal_names = [temp_metal_names0[4], temp_metal_names0[5], temp_metal_names0[1], temp_metal_names0[2], temp_metal_names0[3]]

    for ii in range(len(metallicity_names)):
        err_metal, xpeak_metal = compute_onesig_pdf(metallicity_names[ii], combined_metallicity[ii], usepeak=True, silent=True, verbose = True)
        
        #Filling in dictionaries 
        metal_error[metal_err[ii]] = err_metal
        metal_xpeak[xpeak_key[ii]] = xpeak_metal

        #Editing TM_File to add errors for metallicities
        #err_mT = err_metal.transpose()
        #colu_name_idx = TM_tab.index_column(metal_names[ii])
        #c1 = Column(err_tT[0], name=metal_names[ii]+'_Low_Error')
        #c2 = Column(err_tT[1], name=metal_names[ii]+'_High_Error')
        #TM_tab.add_columns([c1, c2], indexes=[colu_name_idx,colu_name_idx])
        TM_tab[metal_names[ii]][detection] = xpeak_metal 

    asc.write(TM_tab, TM_new_name, format = 'fixed_width_two_line')
    np.savez(path+'Te_propdist_dict.npz', **Te_propdist_dict)   #error from compute one sig
    np.savez(path+'Te_xpeaks.npz', **Te_xpeaks)
    np.savez(path+'metal_errors.npz', **metal_error)
    np.savez(path+'metal_xpeaks.npz', **metal_xpeak)
    np.savez(path+'metallicity_pdf.npz', **metallicity_propdist_dict)
    np.savez(path+'flux_propdist.npz', **flux_propdist_dict)
    np.savez(path+'flux_errors.npz',**flux_lowhigh)
    np.savez(path+'flux_xpeak.npz', **flux_xpeak)
    np.savez(path+'Te_errors.npz',**Te_lowhigh)
    

    
def plotting_linear_errors(path, flux_file):
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
    
