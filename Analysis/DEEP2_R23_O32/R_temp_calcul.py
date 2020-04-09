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
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table


#For generalizing for several users
from getpass import getuser
from astropy import units as u

from Metallicity_Stack_Commons.Metallicity_Stack_Commons.temp_metallicity_calc import \
    R_calculation, temp_calculation, metallicity_calculation

from Metallicity_Stack_Commons.Metallicity_Stack_Commons import fitspath_reagen as fitspath_ini
from Metallicity_Stack_Commons.Metallicity_Stack_Commons import k_dict

k_4363  = k_dict['OIII_4363']
k_5007  = k_dict['OIII_5007']
k_3727  = k_dict['OII_3727']
k_4959  = k_dict['OIII_4958']
k_HBETA = k_dict['HBETA']

#Constants
a = 13205
b = 0.92506
c = 0.98062

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
    ID = combine_fits['bin_ID'].data


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
    if dustatt ==True:
        dust_attenuation(fitspath, combine_flux_ascii)
        print('Using attenuated dust values')
        atten_val = asc.read(dust_ascii)

        EBV = atten_val['E(B-V)']

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


    ###Calls for the stacked measurements for R23_O32 project 

    print('Running R, Temperature, and Metallicity Calculations for Stacked Spectra')
    #Ascii Table from FITTING
    OIII5007     = combine_fits['OIII_5007_Flux_Observed'].data
    OIII4959     = combine_fits['OIII_4958_Flux_Observed'].data
    raw_OIII4363 = combine_fits['OIII_4363_Flux_Observed'].data
    Hgamma       = combine_fits['HGAMMA_Flux_Observed'].data
    HBETA        = combine_fits['HBETA_Flux_Observed'].data
    OII3727      = combine_fits['OII_3727_Flux_Observed'].data
    R23_avg      = combine_fits['logR23_avg'].data
    O32_avg      = combine_fits['logO32_avg'].data
    N_Galaxy     = combine_fits['N_stack'].data
    ID           = combine_fits['bin_ID'].data

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
    R_value= R_calculation(OIII4363, OIII5007, EBV)
    T_e= temp_calculation(R_value)
    metal_dict = metallicity_calculation(T_e, Two_Beta, Three_Beta)


    
    n=  ('bin_ID', 'Detection', 'logR23_Composite', 'logO32_Composite', 'logR23_avg', 'logO32_avg', 'N_stack', 'OIII_5007_Flux_Observed','OIII_5007_S/N', 'OIII_4959_Flux_Observed','OIII_4959_S/N','OIII_4363_Flux_Observed','OIII_4363_S/N', 'HBETA_Flux_Observed','HBETA_S/N', 'OII_3727_Flux_Observed','OII_3727_S/N','T_e', 'O+/H', 'O++/H','12+log(O/H)', 'log(O+/H)', 'log(O++/H)')

    '''c_add = Column(indicate, name='Detection')
    combine_fits.add_columns(c_add) #index = 0 should add the detection column to the front of the combine_flux table'''

    tab0 = Table([ID , indicate, R23_composite, O32_composite, R23_avg, O32_avg, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, T_e, metal_dict['O+/H'], metal_dict['O++/H'], metal_dict['12+log(O/H)'], metal_dict['log(O+/H)'], metal_dict['log(O++/H)']], names=n)
    asc.write(tab0, out_ascii, format='fixed_width_two_line')

    
    tab0.write(out_fits,format = 'fits', overwrite = True)


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



