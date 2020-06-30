import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from astropy.table import Table
import os
from os.path import exists
import glob
from datetime import date



def check_verification_table(fitspath_ini, dataset, combine_flux_ascii):
    verification_table = fitspath_ini+'verification_tables/'+dataset+'_verification_tab.tbl'
    if exists(verification_table):
        return verification_table
    else:
        print('Making verification table')
        verification_tables.verification_master(fitspath,dataset, combine_flux_ascii)
        return verification_table


###TWO FUNCTIONS THAT CREATE THE ACTIVE AND MASTER VERIFICATION TABLES
###RUN EVERYTIME THE GENERAL FUNCTIONS ARE CALLED

# Writes Ascii Files with all the Correct data in them

import numpy as np
import matplotlib.pyplot as plt
# import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ
import glob


def verification_master(fitspath_ini, dataset, combine_flux_ascii):
    # temp_tab = asc.read(temp_table)
    com_flux = asc.read(combine_flux_ascii)

    ver_tab = fitspath_ini + 'verification_tables/bin_validation.tbl'

    ID = com_flux['bin_ID'].data
    raw_OIII4363 = com_flux['OIII_4363_Flux_Observed'].data
    SN_4363 = com_flux['OIII_4363_S/N'].data
    N_Galaxy = com_flux['N_stack'].data
    Hgamma_SN = com_flux['HGAMMA_S/N'].data
    Hgamma = com_flux['HGAMMA_Flux_Observed'].data

    indicate = np.zeros(len(ID))
    OIII_4363 = np.zeros(len(ID))
    up_limit = (Hgamma / Hgamma_SN) * 3

    # Detection Where Statements
    if dataset == 'Voronoi20': det_4363 = np.where((ID == 0) | (ID == 2) | (ID == 3) | (ID == 5) | (ID == 6))[0]
    if dataset == 'Voronoi14': det_4363 = np.where((ID == 0) | (ID == 7) | (ID == 10) | (ID == 11) | (ID == 12))[0]
    if dataset == 'Voronoi10': det_4363 = np.where((ID == 1) | (ID == 9) | (ID == 18) | (ID == 21))[0]
    if dataset == 'Grid': det_4363 = np.where((ID == 11) | (ID == 13) | (ID == 19) | (ID == 20) | (ID == 21))[0]
    if dataset == 'R23_Grid': det_4363 = np.where((ID == 0) | (ID == 4) | (ID == 5) | (ID == 6))[0]
    if dataset == 'O32_Grid': det_4363 = np.where((ID == 6))[0]
    if dataset == 'Double_Bin': det_4363 = \
    np.where((ID == 0) | (ID == 1) | (ID == 2) | (ID == 7) | (ID == 9) | (ID == 10) | (ID == 11) | (ID == 13))[0]
    # if dataset== 'Double_Bin': det_4363 = np.where((SN_4363>=3) | (ID ==0)| (ID ==1)| (ID ==2) | (ID ==7) | (ID ==9) | (ID ==10) | (ID ==11) | (ID ==13))[0]
    if dataset == 'n_Bins':
        det_4363 = np.where((ID == 10) | (ID == 11) | (ID == 14) | (ID == 15) | (ID == 20) | (ID == 23) | (ID == 26))[0]
        rlimit = \
        np.where((ID == 5) | (ID == 7) | (ID == 8) | (ID == 13) | (ID == 16) | (ID == 17) | (ID == 19) | (ID == 22))[0]

    print('det_4363:', det_4363)
    # 5, 7, 8, 10, 11, 13, 14, 15, 16, 17, 19, 20, 22, 23,26
    # 2,7, 8, 10, 14, 15, 20, 23, 26

    indicate[det_4363] = 1
    indicate[rlimit] = 0.5

    OIII_4363[det_4363] = raw_OIII4363[det_4363]
    OIII_4363[rlimit] = up_limit[rlimit]
    nan_det = np.where((indicate == 0))[0]
    OIII_4363[nan_det] = up_limit[nan_det]
    print('indicate', indicate)

    n = ('bin_ID', 'Detection', 'OIII_4363_Flux_Observed', 'OIII_4363_S/N')
    tab1 = Table([ID, indicate, OIII_4363, SN_4363], names=n)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')


###Not Being Used/Notes
###Equivalent Width correction
'''bottom side of the iceburg / continum 
take negative component of gauss and subtract off the postive component 

total gauss - the positve component? 
x = double gaussian of o1[0,1,2] = 0 
x-o1[3] from [-2.5sigma to 2.5 sigma]
equivalent width in terms of angstrumns 
update plots    O_s_ion, O_d_ion,  O_s_ion', 'O_d_ion' '''

'''det_O32 = O32_all[det_4363]
det_R23 = R23_all[det_4363]
det_OH  = com_O_log[det_4363]
det_ID  = ID[det_4363]
det_SN_4363 = SN_4363[det_4363]
det_obs4363 = all_4363[det_4363]
det_N_Gal = N_Gal_all[det_4363]
det_temp = temp_all[det_4363]


n =('ID','Detection', 'R23_Composite', 'O32_Composite', 'R_23_Average', 'O_32_Average', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959', 'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727', 'S/N_3727', 'Temperature', 'com_O_log')
 tab1 = Table([ID, indicate, R23_all, O32_all, R23_avg, O32_avg, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII_4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, temp_all, com_O_log], names=n)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')'''

# OIII5007 = com_flux['OIII_5007_Flux_Observed'].data
# OIII4959 = com_flux['OIII_4958_Flux_Observed'].data

# Hgamma = com_flux['HGAMMA_Flux_Observed'].data
# HBETA    = com_flux['HBETA_Flux_Observed'].data
# OII3727  = com_flux['OII_3727_Flux_Observed'].data
# R23_avg      = com_flux['R_23_Average'].data
# O32_avg      = com_flux['O_32_Average'].data


# SN_Hgamma    = com_flux['HGAMMA_S/N'].data
# SN_5007       = com_flux['OIII_5007_S/N'].data
# SN_4959       = com_flux['OIII_4958_S/N'].data
# SN_4363       = com_flux['OIII_4363_S/N'].data
# SN_HBETA      = com_flux['HBETA_S/N'].data
# SN_3727       = com_flux['OII_3727_S/N'].data


# O32_all = temp_tab['O32_Composite'].data
# R23_all = temp_tab['R23_Composite'].data
# com_O_log = temp_tab['com_O_log']  #This is the 12+log(OH) value

# all_4363 =temp_tab['Observed_Flux_4363'].data
# N_Gal_all = temp_tab['N_Galaxies'].data
# temp_all = temp_tab['Temperature'].data


'''
###Verified Sources###
ver_20 = np.array([])
ver_14 = np.array([0,7,10,11,12])
ver_10 = np.array([1,9,18,21])
ver_Grid = np.array([11,13,19,20,21])
ver_sR23 = np.array([0,4,5,6])
ver_sO32 = np.array([6])'''

'''fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/'
dataset = 'Double Bin'
temp_table = fitspath+ '/'+dataset+'_temperatures_metalicity.tbl'
combine_flux_ascii = fitspath+dataset+'_combined_flux_table.tbl'
temp_tab = asc.read(temp_table)
com_flux = asc.read(combine_flux_ascii)
ver_tab = fitspath+'/'+dataset+'_verification.tbl'''
