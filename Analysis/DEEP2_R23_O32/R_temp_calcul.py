
###THIS CODE IS RUN IN THE GENERAL FUNCTION###

#Calculates the R value, electron temperature, and metallicity from the flux table
#produced by the zoom_and_gauss_general functions


##############Functions#######################
# R_calculation(OIII4363, OIII5007, OIII4959, EBV, k_4363, k_5007)
# temp_calculation(R)
# metallicity_calculation(T_e, TWO_BETA, THREE_BETA)
# limit_function(combine_flux_ascii)
# run_function(fitspath, dataset, out_ascii='', out_fits='', pdf_name='',  combine_flux_ascii='', dust_ascii='', dustatt= False)
'''
Input variables: 
fitspath
dataset
out_ascii          -> name of the output ascii table that file produces
out_fits           -> name of the output fits table that file produces
pdf_name           -> name of the pdf file with the graphs that the function will produce
combine_flux_ascii -> ascii table with all the emission lines 
dust_ascii         -> name of ascii table with dust attenuation values 
dustatt            -> True/False input; if True dust attenuation values are used for calculations; automatic = false 




Calling order: 
verification tables   --> need to work on making these tables; need to put the verification table into the call
dust attenuation      --> called in function by True or False, but need to pass the table into the function
Called DEEP2 and MACT Data
Depending on which combine_fits table is passed in --> run individual or stacked spectra and makes a table
Plots (currently commented out)
'''
# dust_attenuation(combine_ascii)
# call_cardelli(lam0)
# dust_vs_nondust_table
'''
creates a table for the non-dust and dust attenuation metallicity values
called independently from the run_function 
'''

#Currently running: Grid
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
from mpl_toolkits.mplot3d import Axes3D
import sys

#For generalizing for several users
from getpass import getuser
from astropy import units as u

from Metallicity_Stack_Commons.temp_metallicity_calc import \
    R_calculation, temp_calculation, metallicity_calculation

from Metallicity_Stack_Commons import fitspath_reagen as fitspath_ini

#Constants

a = 13205
b = 0.92506
c = 0.98062

def R_calculation(OIII4363, OIII5007, EBV, k_4363, k_5007):
  
    R_value = OIII4363/(OIII5007*(1+1/3.1))* 10**(0.4*EBV*(k_4363-k_5007))
    return R_value  

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)     
    print(T_e)  
    return T_e  


def metallicity_calculation(T_e, TWO_BETA, THREE_BETA):   
    
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17
    
    t_3 = T_e*1e-4
    t_2 = 0.7*t_3 +0.17
    x2 = 1e-4 * 1e3 * t_2**(-0.5)

    O_s_ion_log = np.log10(TWO_BETA) +5.961 +1.676/t_2 - 0.4*np.log10(t_2) - 0.034*t_2 + np.log10(1+1.35*x2)-12
    O_d_ion_log = np.log10(THREE_BETA)+6.200+1.251/t_3 - 0.55*np.log10(t_3) - 0.014*(t_3)-12

    O_s_ion = 10**(O_s_ion_log)
    O_d_ion = 10**(O_d_ion_log)
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) +12

    return O_s_ion , O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log

def limit_function(combine_flux_ascii):
    
    combine_fits= asc.read(combine_flux_ascii)

    Hgamma = combine_fits['HGAMMA_Flux_Observed'].data
    Hgamma_SN    = combine_fits['HGAMMA_S/N'].data

    up_temp = (Hgamma/Hgamma_SN) *3

    #print 'up_temp', up_temp
    return up_temp
    
def run_function(fitspath, dataset, verification_table, out_ascii='', out_fits='', pdf_name='',  combine_flux_ascii='', dust_ascii='', dustatt= False):
    print(combine_flux_ascii)

    ###Combine_Flux_ascii table import 
    combine_fits= asc.read(combine_flux_ascii)
    ID = combine_fits['ID'].data


    #####Verification Table Import#######
    #print('Using verification table' + verification_table)
    ver_tab = asc.read(verification_table)
    ver_detection = ver_tab['Detection']
    ver_detect = np.where((ver_detection ==1))[0]
    ver_rlimit = np.where((ver_detection ==0.5))[0]
    nan_detect = np.where((ver_detection == 0))[0]
    
    detect_ID = ID[ver_detect]
    rlimit_ID = ID[ver_rlimit]

    #Dust Attenuation Values Call
    ###Every time I change the binning method I need to go back and recalculate the dust attentuation##
    if dustatt ==False:
        
        EBV = np.zeros(len(ID))
        k_3727 = np.zeros(len(ID))
        k_HDELTA = np.zeros(len(ID))
        k_Hgamma = np.zeros(len(ID))
        k_HBETA = np.zeros(len(ID))
        k_4363 = np.zeros(len(ID))
        k_4959 = np.zeros(len(ID))
        k_5007 = np.zeros(len(ID))

    if dustatt ==True:
        dust_attenuation(fitspath, combine_flux_ascii)
        print('Using attenuated dust values')
        atten_val = asc.read(dust_ascii)

        EBV = atten_val['E(B-V)']
        k_3727 = atten_val['k_3727']
        k_HDELTA = atten_val['k_HDELTA']
        k_Hgamma = atten_val['k_Hgamma']
        k_HBETA = atten_val['k_HBETA']
        k_4363 = atten_val['k_4363']
        k_4959 = atten_val['k_4958']
        k_5007 = atten_val['k_5007']
    
    #Fits Table Calls
    ###DEEP2 and MACT Data#####
    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')

    #DEEP2 Derived 
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_Te = derived['Te'].data
    der_OH = derived['OH'].data
    ID_der = derived['ID'].data
    #der_OH_log = np.log10(er_OH_log)

    #MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_Te_MACT = derived_MACT['Te'].data
    der_OH_MACT = derived_MACT['OH'].data
    ID_der_MACT = derived_MACT['ID'].data
    


    
    
    ###Calls for individual detections for R23_O32/Reagen project 
    if 'two_beta' in combine_fits.keys():
        print('Running Metallicity Calculation for Individual Spectra')
        OIII5007 = combine_fits['OIII5007'].data
        OIII4959 = combine_fits['OIII4959'].data
        HBETA    = combine_fits['HBeta'].data
        HGamma = combine_fits['HGAMMA'].data
        SNR_HG = combine_fits['SNR_HG'].data
        raw_OIII4363 = combine_fits['OIII4363'].data
        SNR_4363 = combine_fits['SNR_4363'].data
        R23_individual = combine_fits['Individual_R23'].data
        O32_individual = combine_fits['Individual_O32'].data
        Detection = combine_fits['Detection'].data

        up_limit = (Hgamma/SN_Hgamma) *3
        print('up_limit', up_limit)
        
        OIII4363 = np.zeros(len(raw_OIII4363))
        indicate = np.zeros(len(raw_OIII4363))
        for ii in range(len(OIII4363)):
            print(SN_4363[ii])
            if ver_detection[ii] == 1: 
                OIII4363[ii]= raw_OIII4363[ii]
                indicate[ii]= 1
            if ver_detection[ii] == 0.5:
                OIII4363= up_limit
                indicate= 0.5
            if ver_detction[ii] ==0: 
                print('upper limit')
                print('4363: ', up_limit[ii])
                OIII4363[ii]= up_limit[ii]
                indicate[ii]= 0
            print(OIII4363)
            print(indicate) 
        
        Two_Beta = combine_fits['two_beta'].data
        Three_Beta = combine_fits['three_beta'].data
        T_e = combine_fits['Temperature'].data
        Source_ID = combine_fits['Source_ID'].data

        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, Two_Beta, Three_Beta)

        
        n=  ('Source_ID', 'R23', 'O32', 'Observed_Flux_5007', 'Observed_Flux_4959', 'Flux_4363','Observed_Flux_HBETA','Temperature', 'Detection','O_s_ion', 'O_d_ion', 'com_O_log','log_O_s', 'log_O_d' )  #, '3727_HBETA', '5007_HBETA', '4959_HBETA', '5007_3727', '4959_3727', '4363_5007')
        
        tab0 = Table([Source_ID, R23_individual, O32_individual, OIII5007, OIII4959, OIII4363, HBETA, T_e, Detection, O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d], names=n)   # 3727_HBETA,5007_HBETA,4959_HBETA,5007_3727, 4959_3727, 4363_5007








    ####Calls for the individual detections for the MassLuminosity/Caroline project 
    if 'Log10(Mass)' in combine_fits.keys():
        print('Running Metallicity Calculation for Individual Spectra')
        OIII5007 = combine_fits['OIII5007_Flux'].data
        OIII4959 = combine_fits['OIII4959_Flux'].data
        HBETA    = combine_fits['HBETA_Flux'].data
        Mass = combine_fits['Log10(Mass)'].data
        LHbeta = combine_fits['HBeta_Luminosity'].data
        OII = combine_fits['OII_Flux'].data

        #Two_Beta = combine_fits['two_beta'].data
        #Three_Beta = combine_fits['three_beta'].data
        T_e = combine_fits['Te'].data
        Source_ID = combine_fits['OBJNO'].data
        Two_Beta = OII/HBETA
        Three_Beta = (OIII5007+OIII4959)/HBETA

        
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, Two_Beta, Three_Beta)

        
        n=  ('Source_ID', 'Mass', 'LHbeta', 'Observed_Flux_5007', 'Observed_Flux_4959', 'Observed_Flux_HBETA', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log', 'log_O_s', 'log_O_d')  #, '3727_HBETA', '5007_HBETA', '4959_HBETA', '5007_3727', '4959_3727', '4363_5007')
        
        tab0 = Table([Source_ID, Mass, LHbeta, OIII5007, OIII4959, HBETA, T_e, O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d], names=n)   # 3727_HBETA,5007_HBETA,4959_HBETA,5007_3727, 4959_3727, 4363_5007 






    ###Calls for the stacked measurements for R23_O32 project 
    else:
        print('Running R, Temperature, and Metallicity Calculations for Stacked Spectra')
        #Ascii Table from FITTING 
        OIII5007     = combine_fits['OIII_5007_Flux_Observed'].data
        OIII4959     = combine_fits['OIII_4958_Flux_Observed'].data
        raw_OIII4363 = combine_fits['OIII_4363_Flux_Observed'].data
        Hgamma       = combine_fits['HGAMMA_Flux_Observed'].data
        HBETA        = combine_fits['HBETA_Flux_Observed'].data
        OII3727      = combine_fits['OII_3727_Flux_Observed'].data
        R23_avg      = combine_fits['R_23_Average'].data
        O32_avg      = combine_fits['O_32_Average'].data
        N_Galaxy     = combine_fits['N_Galaxies'].data
        ID           = combine_fits['ID'].data

        SN_Hgamma    = combine_fits['HGAMMA_S/N'].data
        SN_5007      = combine_fits['OIII_5007_S/N'].data
        SN_4959      = combine_fits['OIII_4958_S/N'].data
        SN_4363      = combine_fits['OIII_4363_S/N'].data
        SN_HBETA     = combine_fits['HBETA_S/N'].data
        SN_3727      = combine_fits['OII_3727_S/N'].data

        R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
        O32_composite = np.log10((1.33*OIII5007)/OII3727)

        print('R23_composite', R23_composite)
        print('O32_composite', O32_composite)
    
        up_limit = (Hgamma/SN_Hgamma) *3
        print('up_limit', up_limit)
        
        OIII4363 = np.zeros(len(raw_OIII4363))
        indicate = np.zeros(len(raw_OIII4363))

        for ii in range(len(OIII4363)):
            
            if ver_detection[ii] == 1:
                OIII4363[ii]= raw_OIII4363[ii]
                indicate[ii]= 1 
            if ver_detection[ii] == 0.5:
                OIII4363[ii]= up_limit[ii] 
                indicate[ii]= 0.5
            if ver_detection[ii] ==0: 
                OIII4363[ii]= up_limit[ii]
                
        

        #Line Ratios
        O3727_HBETA = OII3727/HBETA
        O5007_HBETA = OIII5007/HBETA
        O4959_HBETA = OIII4959/HBETA
        O4363_O5007 = OIII4363/OIII5007
        O4363_O4959 = OIII4363/OIII4959
    
        O5007_O3727 = OIII5007/OII3727
        O4959_O3727 = OIII4959/OII3727
        
        #Attenuated Ratios
        der_4363_5007  = O4363_O5007 * 10**(0.4*EBV*(k_4363-k_5007))
        der_4363_4959  = O4363_O4959 * 10**(0.4*EBV*(k_4363-k_4959))
        der_3727_HBETA = O3727_HBETA * 10**(0.4*EBV*(k_3727-k_HBETA))
        der_4959_HBETA = O4959_HBETA * 10**(0.4*EBV*(k_4959-k_HBETA))
        der_5007_HBETA = O5007_HBETA * 10**(0.4*EBV*(k_5007-k_HBETA))

        if dustatt == False: 
            Two_Beta = OII3727/HBETA
            Three_Beta= (OIII5007* (1+1/3.1))/HBETA
        else:
            Two_Beta = der_3727_HBETA                    
            Three_Beta= der_5007_HBETA*(1+1/3.1) 
    
        #Raw Data
        R_value= R_calculation(OIII4363, OIII5007, EBV, k_4363, k_5007)  
        T_e= temp_calculation(R_value)  
        O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metallicity_calculation(T_e, Two_Beta, Three_Beta)



        n=  ('ID', 'Detection', 'R23_Composite', 'O32_Composite', 'R_23_Average', 'O_32_Average', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959', 'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727', 'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log', 'log_O_s','log_O_d')
        #, '3727_HBETA', '5007_HBETA', '4959_HBETA', '5007_3727', '4959_3727', '4363_5007')

         
        '''c_add = Column(indicate, name='Detection')
        combine_fits.add_columns(c_add) #index = 0 should add the detection column to the front of the combine_flux table'''
        



        tab0 = Table([ID , indicate, R23_composite, O32_composite, R23_avg, O32_avg, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, T_e, O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d], names=n)   


    
        

    

    asc.write(tab0, out_ascii, format='fixed_width_two_line')

    
    tab0.write(out_fits,format = 'fits', overwrite = True)

    #Plots
    
    ###Getting the the upper/lower limits for the metallicity and temperatures###
    metallicity_npz = np.load(fitspath_ini +'Error_propagation/metal_errors.npz')
    npz_comOlog = metallicity_npz['com_O_log_lowhigh_error']
    metal_error = np.transpose(npz_comOlog)
    print("metal_error:", metal_error)

    npz_te_err = np.load(fitspath_ini + 'Error_propagation/Te_errors.npz')
    npz_te = npz_te_err['Te_lowhigh_error']
    Te_error = np.transpose(npz_te)


    ###Starting the plots 
    pdf_pages = PdfPages(fitspath+pdf_name)
    print(fitspath+pdf_name)
    color_len = len(R23_composite)
    
    color_arr = plt.cm.get_cmap('rainbow')
    
    
    ###Defining the logarithms of values
    lTe = np.log10(T_e)
    lder_Te = np.log10(der_Te)
    lMACT_Te = np.log10(der_Te_MACT)
    l_Te_err = np.log10(Te_error)



    #################################Temperature vs R23_composite #################################
    '''fig1, ax1 = plt.subplots()
    ax1.scatter(T_e[ver_detect], R23_composite[ver_detect], marker = '.', s = 50, color = 'b')
    #ax1.scatter(lTe[nan_detect], R23_composite[nan_detect], marker = '<', s = 25, color = 'b')
    ax1.scatter(T_e[ver_rlimit], R23_composite[ver_rlimit], marker = '<', s = 35, color = 'b')
    
    ax1.errorbar(lTe[ver_detect], R23_composite[ver_detect], xerr = l_Te_err[:,:], ls = 'none')
    for aa in range(len(ID)): ax1.annotate(ID[aa], (lTe[aa], R23_composite[aa]), fontsize = '6')

    ax1.scatter(lder_Te, der_R23, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for bb in range(len(ID_der)): ax1.annotate(ID_der[bb], (lder_Te[bb], der_R23[bb]), fontsize ='2')

    ax1.scatter(lMACT_Te, der_R23_MACT, s = 20, marker = 'P', color = 'r', edgecolors = 'None', alpha = 0.5)
    for qq in range(len(ID_der_MACT)): ax1.annotate(ID_der_MACT[qq], (lMACT_Te[qq], der_R23_MACT[qq]), fontsize= '2')

    ax1.set_xlabel('Log Temperature (K)')
    ax1.set_ylabel('R_23')
    ax1.set_title('Temperatures_vs_R23')
    fig1.savefig(pdf_pages, format = "pdf")
    '''
    
###########################################################################################
##############################Temperature vs R23_composite With Limits#####################
    fig5, ax5 = plt.subplots()
    
    ax5.scatter(T_e[ver_detect], R23_composite[ver_detect], marker = '.',s = 50, color = 'b')
    #ax5.scatter(T_e[nan_detect], R23_composite[nan_detect], marker = '<', s = 25, color = 'b')
    ax5.scatter(T_e[ver_rlimit], R23_composite[ver_rlimit], marker = '<',s = 35, color = 'b')
    
    ax5.errorbar(T_e[ver_detect], R23_composite[ver_detect], xerr = Te_error[:,:], ls = 'none')
    for xx in ver_detect:ax5.annotate(ID[xx], (T_e[xx], R23_composite[xx]), fontsize = '6')
    for xx in ver_rlimit:ax5.annotate(ID[xx], (T_e[xx], R23_composite[xx]), fontsize = '6')

    ax5.scatter(der_Te, der_R23, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for b in range(len(ID_der)): ax5.annotate(ID_der[b], (der_Te[b], der_R23[b]), fontsize = '2')

    ax5.scatter(der_Te_MACT, der_R23_MACT, s =20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for q in range(len(ID_der_MACT)): ax5.annotate(ID_der_MACT[q], (der_Te_MACT[q], der_R23_MACT[q]), fontsize= '2')

    ax5.set_xlabel('Temperature (K)')
    ax5.set_ylabel('R_23')
    ax5.set_title('Temperatures_vs_R23 with Limits on Temperature')
    ax5.set_xlim(5000,21500)
    fig5.savefig(pdf_pages, format = "pdf")

############################################################################################
###########################Temperature vs O32_composite#####################################
    '''fig2, ax2 = plt.subplots()
    
    ax2.scatter(T_e[ver_detect], O32_composite[ver_detect], marker = '.',s = 50, color = 'b')
    #ax2.scatter(lTe[nan_detect], O32_composite[nan_detect], marker = '<',s = 25, color = 'b')
    ax2.scatter(T_e[ver_rlimit], O32_composite[ver_rlimit], marker = '<',s = 35, color = 'b')
    
    ax2.errorbar(lTe[ver_detect], O32_composite[ver_detect], xerr = l_Te_err[:,:], ls = 'none')
    for cc in range(len(ID)): ax2.annotate(ID[cc], (lTe[cc], O32_composite[cc]), fontsize = '6')

    ax2.scatter(lder_Te, der_O32, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for ff in range(len(ID_der)):ax2.annotate(ID_der[ff], (lder_Te[ff], der_O32[ff]), fontsize = '2')

    ax2.scatter(lMACT_Te, der_O32_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5,edgecolors = 'None')
    for ss in range(len(ID_der_MACT)): ax2.annotate(ID_der_MACT[ss], (lMACT_Te[ss], der_O32_MACT[ss]), fontsize= '2')
    
    ax2.set_xlabel('Log Temperature (K)')
    ax2.set_ylabel('O_32')
    ax2.set_title('Temperatures_vs_O32')
    fig2.savefig(pdf_pages, format = 'pdf')
    '''

###############################################################################################
####################################Temperature vs O32_Composite with limits##################
    #With Limits
    fig6, ax6 = plt.subplots()
    
    ax6.scatter(T_e[ver_detect], O32_composite[ver_detect], marker = '.',s=50, color = 'b')
    #ax6.scatter(T_e[nan_detect], O32_composite[nan_detect], marker = '<',s=25, color = 'b')
    ax6.scatter(T_e[ver_rlimit], O32_composite[ver_rlimit], marker = '<',s=35, color = 'b')
    ax6.errorbar(T_e[ver_detect], O32_composite[ver_detect], xerr = Te_error[:,:], ls = 'none')
    for c in ver_detect:ax6.annotate(ID[c], (T_e[c], O32_composite[c]), fontsize = '6')
    for c in ver_rlimit:ax6.annotate(ID[c], (T_e[c], O32_composite[c]), fontsize = '6')

    ax6.scatter(der_Te, der_O32, s=20, marker = '*', color = 'k',edgecolors = 'None')
    for f in range(len(ID_der)):ax6.annotate(ID_der[f], (der_Te[f], der_O32[f]), fontsize = '2')

    ax6.scatter(der_Te_MACT, der_O32_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors ='None')
    for s in range(len(ID_der_MACT)):ax6.annotate(ID_der_MACT[s], (der_Te_MACT[s], der_O32_MACT[s]), fontsize= '2')
    
    ax6.set_xlabel('Temperature (K)')
    ax6.set_ylabel('O_32')
    ax6.set_title('Temperatures_vs_O32 with Limits on Temperature')
    ax6.set_xlim(5000,21500)
    fig6.savefig(pdf_pages, format = 'pdf')
    

################################################################################################
#########################################R23_composite vs metallicity###########################
    fig3, ax3 = plt.subplots()
    
    ax3.scatter(R23_composite[ver_detect], com_O_log[ver_detect], marker = '.', s = 50, color = 'b')
    #ax3.scatter(R23_composite[nan_detect], com_O_log[nan_detect], marker = '^', s = 20, color = 'b')
    ax3.scatter(R23_composite[ver_rlimit], com_O_log[ver_rlimit], marker = '^', s = 35, color = 'b')
    
    ax3.errorbar(R23_composite[ver_detect], com_O_log[ver_detect], yerr = metal_error[:,:], ls = 'none')
    for zz in ver_detect:ax3.annotate(ID[zz], (R23_composite[zz],com_O_log[zz]), fontsize = '6')
    for zz in ver_rlimit:ax3.annotate(ID[zz], (R23_composite[zz],com_O_log[zz]), fontsize = '6')
    
    ax3.scatter(der_R23, der_OH, s= 20, marker = '*', color = 'k', edgecolors = 'None')
    for gg in range(len(ID_der)): ax3.annotate(ID_der[gg], (der_R23[gg], der_OH[gg]), fontsize = '2')
    ax3.scatter(der_R23_MACT, der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors='None')
    for g in range(len(ID_der_MACT)): ax3.annotate(ID_der_MACT[g], (der_R23_MACT[g], der_OH_MACT[g]), fontsize= '2')

    ax3.set_xlabel('R23')
    ax3.set_ylabel('12+log(O/H) Te')
    ax3.set_title('R23 vs. Composite Metallicity')
    fig3.savefig(pdf_pages, format = 'pdf')

##################################################################################################
########################################O32_composite vs Metallicity##############################
    fig4, ax4 = plt.subplots()
    
    ax4.scatter(O32_composite[ver_detect], com_O_log[ver_detect], marker = '.',s =50, color = 'b')
    #ax4.scatter(O32_composite[nan_detect], com_O_log[nan_detect], marker = '^',s =20, color = 'b')
    ax4.scatter(O32_composite[ver_rlimit], com_O_log[ver_rlimit], marker = '^',s =35, color = 'b')
    
    ax4.errorbar(O32_composite[ver_detect], com_O_log[ver_detect], yerr = metal_error[:,:], ls = 'none')
    for ww in ver_detect:ax4.annotate(ID[ww], (O32_composite[ww], com_O_log[ww]), fontsize = '6')
    for ww in ver_rlimit:ax4.annotate(ID[ww], (O32_composite[ww], com_O_log[ww]), fontsize = '6')
    
    ax4.scatter(der_O32,der_OH, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for hh in range(len(ID_der)): ax4.annotate(ID_der[hh], (der_O32[hh], der_OH[hh]), fontsize = '2')

    ax4.scatter(der_O32_MACT,der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for h in range(len(ID_der_MACT)): ax4.annotate(ID_der_MACT[h], (der_O32_MACT[h], der_OH_MACT[h]), fontsize= '2')

    ax4.set_xlabel('O32')
    ax4.set_ylabel('12+log(O/H) Te')
    ax4.set_title('O32 vs. Composite Metallicity')
    
    fig4.savefig(pdf_pages, format ='pdf')


    

    fig7, ax7 = plt.subplots()
    x_value = np.log10(OIII4363/OIII5007)
    y_value = np.log10(T_e)
    z_value = np.log10(der_Te)
    
    
    ax7.scatter(x_value[ver_detect], y_value[ver_detect], marker = 'o', s= 30,  color = 'b')
    #ax7.scatter(x_value[nan_detect], y_value[nan_detect], marker = 'x', s= 30,  color = 'b')
    ax7.scatter(x_value[ver_rlimit], y_value[ver_rlimit], marker = 'o', s= 30,  color = 'b')
    for r in ver_detect: ax7.annotate(ID[r], (x_value[r], y_value[r]))
    for r in ver_rlimit: ax7.annotate(ID[r], (x_value[r], y_value[r]))
    xxx = np.arange(min(x_value),max(x_value), 1)
    xx_value = np.zeros(len(xxx))
    
    t_E_value = np.zeros(len(xxx))
    for xx in range(len(xxx)):
        xx_value[xx] = xxx[xx]
        T_e =  a*(-np.log10(xxx[xx])-b)**(-1*c)
        t_E_value[xx] = T_e
    ax7.plot(xx_value, t_E_value, 'k')
    
    ax7.set_xlabel('log(OIII4363/OIII5007)')
    ax7.set_ylabel('log(T_e)')
    ax7.set_title('log(OIII4363/OIII5007) vs log(T_e)')
    fig7.savefig(pdf_pages, format = 'pdf')

    pdf_pages.close()

    
'''
    #3D plots
    fig_3d= plt.figure(figsize=(10,8))
    ax_3d = plt.axes(projection='3d')
    ax_3d.set_xlabel('R_23')
    ax_3d.set_ylabel('O_32')
    ax_3d.set_zlabel('Temperature')
    ax_3d.set_zlim(4000,26000)
    ax_3d.scatter(R23_composite, O32_composite, T_e, marker='.', linewidths = None)
    #plt.show()


    #print 'indicate:', indicate 
    #print 'mDect:', mDect
    Te_marker = []
    #print Te_marker
    #print('!!!!!', len(OIII4363))
    for oo in range(len(OIII4363)):
        #print indicate[oo]
        #print mDect[oo]
        if indicate[oo]== 0 and mDect[oo]== 0: Marker = '<'
        if indicate[oo]== 1 and mDect[oo]== 0: Marker = '<'
        if indicate[oo]== 0 and mDect[oo]== 1: Marker = '.'
        if indicate[oo]== 1 and mDect[oo]== 1: Marker = '.'
        Te_marker.append(Marker)
    #print Te_marker

    M_marker = []
    #print M_marker
    for ww in range(len(OIII4363)):
        #print indicate[ww]
        #print mDect[ww]
        if indicate[ww]== 0 and mDect[ww]== 0: Marker= '^'
        if indicate[ww]== 1 and mDect[ww]== 0: Marker= '^'
        if indicate[ww]== 0 and mDect[ww]== 1:  Marker = '.'
        if indicate[ww]== 1 and mDect[ww]== 1:  Marker = '.'
        M_marker.append(Marker)
    #print M_marker









    fig1.clear()
    fig2.clear()
    fig3.clear()
    fig4.clear()
    fig5.clear()
    fig6.clear()
    fig7.clear()
    #fig_3d.clear()
    '''


def dust_attenuation(fitspath, combine_ascii):
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    
    combine_asc = asc.read(combine_ascii)
    ini_con = 0.468
    ID = combine_asc['ID']
    HBeta= combine_asc['HBETA_Flux_Observed'].data
    HGamma= combine_asc['HGAMMA_Flux_Observed'].data
    
    lam0_OII = combine_asc['OII_3727_X_bar'].data
    lam0_HDELTA = combine_asc['HDELTA_X_bar'].data
    lam0_Hgamma = combine_asc['HGAMMA_X_bar'].data
    lam0_HBETA = combine_asc['HBETA_X_bar'].data
    lam0_4363 = combine_asc['OIII_4363_X_bar'].data
    lam0_4958 = combine_asc['OIII_4958_X_bar'].data
    lam0_5007 = combine_asc['OIII_5007_X_bar'].data

    k_3727 = call_cardelli(lam0_OII)
    k_HDELTA = call_cardelli(lam0_HDELTA)
    k_Hgamma = call_cardelli(lam0_Hgamma)
    k_HBETA = call_cardelli(lam0_HBETA)
    k_4363 = call_cardelli(lam0_4363)
    k_4958 = call_cardelli(lam0_4958)
    k_5007 = call_cardelli(lam0_5007)

    
    
    EBV= np.log10((HBeta/HGamma)*(ini_con))*2.5*(1/(k_Hgamma-k_HBETA))
    for nn in range(len(HGamma)):
        if EBV[nn] <= 0: EBV[nn] = 0
    
    #print EBV
    A_3727 = EBV*k_3727
    A_HDELTA = EBV*k_HDELTA
    A_Hgamma = EBV*k_Hgamma
    A_HBETA = EBV*k_HBETA
    A_4363 = EBV*k_4363
    A_4958 = EBV*k_4958
    A_5007 = EBV*k_5007
    #print "A_3727:", A_3727

    out_ascii = fitspath+'/dust_attentuation_values.tbl'
    #if not exists(out_ascii_single):
    n2= ('ID','k_3727', 'k_HDELTA', 'k_Hgamma', 'k_HBETA', 'k_4363', 'k_4958', 'k_5007', 'E(B-V)')
    tab1 = Table([ID,k_3727, k_HDELTA, k_Hgamma, k_HBETA , k_4363, k_4958, k_5007, EBV], names=n2)
    asc.write(tab1, out_ascii, format='fixed_width_two_line')
    
    
def call_cardelli(lam0): #, extrapolate=False):
    #lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]* u.angstrom
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    lambda0 = lam0*u.angstrom
    k_values= cardelli(lambda0,R=3.1)
    return k_values


'''k_ascii = fitspath_ini+'/cardelli_k_values.tbl'
    n3 = ('Line_Name','K_value')
    tab3 = Table([line_name, k_values])
    asc.write(tab3, k_ascii, format = 'fixed_width_two_line')'''




def dust_vs_nondust_table(fitspath, dust_metal_table, nondust_metal_table, dust_atten_values, name):
    #dust_vs_nondust_table = fitspath + 'dust_and_nondust_metaltab.tbl'
    dust_vs_nondust_table = fitspath + name
    
    #Non Dust attentuation
    nondust = asc.read(nondust_metal_table)
    dust = asc.read(dust_atten_values)
    
    nondust_metal = nondust['com_O_log'].data
    dust_metal = dust['com_O_log'].data
    ID = nondust['ID'].data
    R23_composite = nondust['R23_Composite'].data
    O32_composite = nondust['O32_Composite'].data
    Temperature = nondust['Temperature'].data

    n_dust = ('ID','R23_Composite','O32_Composite', 'Non-Dust Attenuated Metallicities','Dust Attenuated Metallicities','Temperature')
    tab_dust = Table([ID, R23_composite, O32_composite, nondust_metal, dust_metal, Temperature], names =n_dust)
    asc.write(tab_dust, dust_vs_nondust_table, format = 'fixed_width_two_line')



