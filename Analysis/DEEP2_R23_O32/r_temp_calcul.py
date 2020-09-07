# THIS CODE IS RUN IN THE GENERAL FUNCTION
# Calculates the R value, electron temperature, and metallicity from the flux table
# produced by the zoom_and_gauss_general functions

import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table
from os.path import join

from Metallicity_Stack_Commons.analysis.temp_metallicity_calc import \
    R_calculation, temp_calculation, metallicity_calculation

from Metallicity_Stack_Commons import fitspath_reagen as fitspath_ini
from Metallicity_Stack_Commons import k_dict
from Metallicity_Stack_Commons.column_names import filename_dict

from Zcalbase_gal.analysis.deep2_r23_o32 import name_dict

k_4363 = k_dict['OIII_4363']
k_5007 = k_dict['OIII_5007']
k_3727 = k_dict['OII_3727']
k_4959 = k_dict['OIII_4958']
k_HBETA = k_dict['HBETA']

# Constants
a = 13205
b = 0.92506
c = 0.98062


def limit_function(combine_flux_ascii):
    """
    Purpose
    Function used to calculate limit of OIII_4363 emission line if no detection
    August 2020: Function no longer used in this file and moved to verification table

    Parameters
    combine_flux_ascii -> table with all the data from emission lines
    """
    combine_fits = asc.read(combine_flux_ascii)

    Hgamma = combine_fits['HGAMMA_Flux_Observed'].data
    Hgamma_SN = combine_fits['HGAMMA_S/N'].data

    up_temp = (Hgamma/Hgamma_SN) * 3

    return up_temp


def run_function(fitspath, dataset, verification_table, dustatt=False):
    """
    Purpose:
    Organize data to calculate R, temperature, and metallicity using MSC and save data into ascii tables
    Input variables:
    fitspath -> location of where the outputted pdf_file will be saved
    dataset  -> keyword used to define binning method
    verification_table -> table used to finalize if bin has detection or not
    dustatt  -> True/False input; if True dust attenuation values are used for calculations; automatic = false

    Calling order:
    verification tables   --> need to work on making these tables; need to put the verification table into the call
    dust attenuation      --> called in function by True or False, but need to pass the table into the function
    Called DEEP2 and MACT Data
    Depending on which combine_fits table is passed in --> run individual or stacked spectra and makes a table

    """
    combine_flux_ascii = join(fitspath, filename_dict['bin_fit'])
    temp_metal_ascii = join(fitspath, filename_dict['bin_derived_prop'])
    temp_metal_revised = join(fitspath, filename_dict['bin_derived_prop_rev'])
    temp_m_gpdf_name = join(dataset, name_dict['temp_metallicity_pdf'])

    # Combine_Flux_ascii table import
    combine_fits = asc.read(combine_flux_ascii)
    id = combine_fits['bin_ID'].data

    # Dust Attenuation
    if dustatt:
        non_atten_value_table = asc.read(temp_metal_ascii)
        EBV = non_atten_value_table['EBV_HgHb']
        out_ascii = join(fitspath, filename_dict['bin_derived_prop_rev_dust'])  # filename_dict['bin_derived_prop_rev']
    else:
        EBV = np.zeros(len(id))
        out_ascii = join(fitspath, filename_dict['bin_derived_prop'])

    # Verification Table Import
    ver_tab = asc.read(verification_table)
    ver_detection = ver_tab['Detection']
    ver_detect = np.where((ver_detection == 1))[0]
    ver_rlimit = np.where((ver_detection == 0.5))[0]

    # Fits Table Calls
    # DEEP2 and MACT Data
    derived = asc.read(fitspath_ini + 'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini + 'MACT_R23_O32_derived.tbl')

    # DEEP2 Derived
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_Te = derived['Te'].data
    der_OH = derived['OH'].data
    ID_der = derived['ID'].data
    # der_OH_log = np.log10(er_OH_log)

    # MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_Te_MACT = derived_MACT['Te'].data
    der_OH_MACT = derived_MACT['OH'].data
    ID_der_MACT = derived_MACT['ID'].data

    # Calls for the stacked measurements for R23_O32 project
    print('Running R, Temperature, and Metallicity Calculations for Stacked Spectra')
    # Ascii Table from FITTING
    OIII5007 = combine_fits['OIII_5007_Flux_Observed'].data
    OIII4959 = combine_fits['OIII_4958_Flux_Observed'].data
    raw_OIII4363 = combine_fits['OIII_4363_Flux_Observed'].data
    Hgamma = combine_fits['HGAMMA_Flux_Observed'].data
    HBETA = combine_fits['HBETA_Flux_Observed'].data
    OII3727 = combine_fits['OII_3727_Flux_Observed'].data
    R23_avg = combine_fits['logR23_avg'].data
    O32_avg = combine_fits['logO32_avg'].data
    N_Galaxy = combine_fits['N_stack'].data
    id = combine_fits['bin_ID'].data

    SN_Hgamma = combine_fits['HGAMMA_S/N'].data
    SN_5007 = combine_fits['OIII_5007_S/N'].data
    SN_4959 = combine_fits['OIII_4958_S/N'].data
    SN_4363 = combine_fits['OIII_4363_S/N'].data
    SN_HBETA = combine_fits['HBETA_S/N'].data
    SN_3727 = combine_fits['OII_3727_S/N'].data

    R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
    O32_composite = np.log10((1.33*OIII5007)/OII3727)

    print('R23_composite', R23_composite)
    print('O32_composite', O32_composite)
    
    up_limit = (Hgamma/SN_Hgamma) * 3
    print('up_limit', up_limit)

    OIII4363 = np.zeros(len(raw_OIII4363))
    indicate = np.zeros(len(raw_OIII4363))

    for ii in range(len(OIII4363)):
        if ver_detection[ii] == 1:
            OIII4363[ii] = raw_OIII4363[ii]
            indicate[ii] = 1
        if ver_detection[ii] == 0.5:
            OIII4363[ii] = up_limit[ii]
            indicate[ii] = 0.5
        if ver_detection[ii] == 0:
            OIII4363[ii] = up_limit[ii]

    # Line Ratios
    O3727_HBETA = OII3727/HBETA
    O5007_HBETA = OIII5007/HBETA
    O4959_HBETA = OIII4959/HBETA
    O4363_O5007 = OIII4363/OIII5007
    O4363_O4959 = OIII4363/OIII4959
        
    # Attenuated Ratios
    der_4363_5007 = O4363_O5007 * 10**(0.4*EBV*(k_4363-k_5007))
    der_4363_4959 = O4363_O4959 * 10**(0.4*EBV*(k_4363-k_4959))
    der_3727_HBETA = O3727_HBETA * 10**(0.4*EBV*(k_3727-k_HBETA))
    der_4959_HBETA = O4959_HBETA * 10**(0.4*EBV*(k_4959-k_HBETA))
    der_5007_HBETA = O5007_HBETA * 10**(0.4*EBV*(k_5007-k_HBETA))

    if dustatt:
        Two_Beta = der_3727_HBETA
        Three_Beta = der_5007_HBETA * (1 + 1 / 3.1)
    else:
        Two_Beta = OII3727 / HBETA
        Three_Beta = (OIII5007 * (1 + 1 / 3.1)) / HBETA

    # Raw Data
    R_value = R_calculation(OIII4363, OIII5007, EBV)
    T_e = temp_calculation(R_value)
    metal_dict = metallicity_calculation(T_e, Two_Beta, Three_Beta)

    n = ('bin_ID', 'Detection', 'logR23_Composite', 'logO32_Composite', 'logR23_avg', 'logO32_avg', 'N_stack',
         'OIII_5007_Flux_Observed', 'OIII_5007_S/N', 'OIII_4959_Flux_Observed', 'OIII_4959_S/N',
         'OIII_4363_Flux_Observed', 'OIII_4363_S/N', 'HBETA_Flux_Observed', 'HBETA_S/N',
         'OII_3727_Flux_Observed', 'OII_3727_S/N', 'T_e', 'O+/H', 'O++/H', '12+log(O/H)', 'log(O+/H)', 'log(O++/H)')
    tab0 = Table([id, indicate, R23_composite, O32_composite, R23_avg, O32_avg, N_Galaxy,
                  OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA,
                  OII3727, SN_3727, T_e, metal_dict['O+/H'], metal_dict['O++/H'], metal_dict['12+log(O/H)'],
                  metal_dict['log(O+/H)'], metal_dict['log(O++/H)']], names=n)
    asc.write(tab0, out_ascii, format='fixed_width_two_line')


"""
This will be used when continuing writing on the paper
    n1=  ('bin_ID', 'Detection', 'R value', 'Electron Temperature', 'O2/HBETA', 'O3/HBETA',
     'O+/H', 'O++/H','12+log(O/H)', 'log(O+/H)', 'log(O++/H)')
    variable_formats=  {'bin_ID': '%i','Detection': '%.1f','R value':'%.3f', 'Electron Temperature': '%.3f', 
    'O2/HBETA': '%.3f', 'O3/HBETA':'%.3f', 'O+/H': '{:.3e}', 'O++/H': '{:.3e}','12+log(O/H)': '%.3f', 
    'log(O+/H)': '%.3f', 'log(O++/H)': '%.3f'}
    tab1 = Table([id, indicate, R_value, Two_Beta, Three_Beta,  T_e, metal_dict['O+/H'], metal_dict['O++/H'],
                  metal_dict['12+log(O/H)'],metal_dict['log(O+/H)'], metal_dict['log(O++/H)']], names=n1)
    asc.write(tab1, '/Users/reagenleimbach/Desktop/Zcalbase_gal/Honors_Thesis/metallicity_table.tex', 
              format='latex', formats= variable_formats)
"""
