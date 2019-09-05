### Finds the temperature and metallicity for individual spectra in a given bin


import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ
import glob
from datetime import date

from . import general

a = 13205
b = 0.92506
c = 0.98062

'''fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/'
dataset = 'n_Bins'
bin_id = 0'''

def ind_detection(fitspath, dataset, bin_id):

    get_det3_tab = asc.read(fitspath+'get_det3_table2.tbl')
    bin_tab = asc.read(fitspath+dataset+'_2d_binning_datadet3.tbl')
    N_gal_tab = asc.read(fitspath+dataset+'_Average_R23_O32_Values.tbl')
    stackmeas_tab = asc.read(fitspath+dataset+'_temperatures_metalicity.tbl')

    #From tables
    Bin_number = bin_tab['Bin_number']
    O2 = get_det3_tab['O2']
    O3 = get_det3_tab['O3']
    Hb = get_det3_tab['Hb']
    N_Galaxies = N_gal_tab['N_Galaxies']
    temp_bin = stackmeas_tab['Temperature']
    ID = N_gal_tab['ID']
    R23 = N_gal_tab['R_23_Average']
    O32 = N_gal_tab['O_32_Average']
    
    #Initializing Arrays
    '''two_beta = np.zeros(int(N_Galaxies[bin_id]))
    three_beta = np.zeros(int(N_Galaxies[bin_id]))
    average_temp = np.zeros(int(N_Galaxies[bin_id]))'''
    
    two_beta = []
    three_beta = []
    average_temp = []
    
    for ii in range(len(O2)):
        if Bin_number[ii] == bin_id:
            print 'Bin_number:', Bin_number[ii], 'O2:', O2[ii], 'O3:', O3[ii], 'Hb:', Hb[ii]
            two_beta.append(O2[ii]/Hb[ii])
            three_beta.append(O3[ii]/Hb[ii])
            average_temp.append(temp_bin[bin_id])

    
    out_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(bin_id)+'_individual_ratios_temp.tbl'
    n = ('two_beta', 'three_beta', 'Temperature')   #'ID', 'R23_Average', 'O32_Average'
    ind_tab = Table([two_beta, three_beta, average_temp], names=n) #ID, R23, O32,
    asc.write(ind_tab, out_ascii, format = 'fixed_width_two_line')


def vstack(fitspath,new_name):
    table_files = glob.glob(fitspath+dataset+'_'+bin_id+'_individual_ratios_temp.tbl')
    table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,new_name, format='fixed_width_two_line', overwrite = True)
    





def ind_metalicity_calculation(T_e,der_3727_HBETA, der_4959_HBETA, der_5007_HBETA, OIII5007, OIII4959, OIII4363, HBETA, OII3727, dustatt = False):
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17

    if dustatt == False: 
        two_beta = OII3727/HBETA
        three_beta= (OIII4959+OIII5007)/HBETA
    else:
        two_beta = der_3727_HBETA                    
        three_beta= der_4959_HBETA+ der_5007_HBETA 
    t_3 = T_e*1e-4
    t_2 = 0.7*t_3 +0.17
    x2 = 1e-4 * 1e3 * t_2**(-0.5)

    O_s_ion_log = np.log10(two_beta) +5.961 +1.676/t_2 - 0.4*np.log10(t_2) - 0.034*t_2 + np.log10(1+1.35*x2)-12
    O_d_ion_log = np.log10(three_beta)+6.200+1.251/t_3 - 0.55*np.log10(t_3) - 0.014*(t_3)-12

    O_s_ion = 10**(O_s_ion_log)
    O_d_ion = 10**(O_d_ion_log)
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) +12

    return O_s_ion , O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log
