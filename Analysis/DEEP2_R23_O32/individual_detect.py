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

#import Metallicity_Stack_Commons
from Metallicity_Stack_Commons import
#from Metallicity_Stack_Commons.analysis.composite_indv_detect import main

a = 13205
b = 0.92506
c = 0.98062

'''fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/'
dataset = 'n_Bins'
bin_id = 0
individual_ascii =  '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(bin_id)+'_individual_ratios_temp.tbl'''


def run_ind_detection(fitspath, dataset, average_value_ascii):
    '''This function runs the function below to create a table with all the galaxies that can be used for individual detections. Each bin then has a table with the individual galaxies/measurements that then gets stacking into one large table. This table can then be passed into R_temp_cal.py to get the individual metallicities and temperatures.

    '''

    N_gal_tab = asc.read(average_value_ascii)  #fitspath+dataset+'_Average_R23_O32_Values.tbl'
    ID = N_gal_tab['bin_ID']
    for aa in range(len(ID)):
        print ID[aa]
        ind_detection(fitspath,dataset,ID[aa])
    new_name = fitspath+'Individual_ratio_temperature.tbl'
    individual_galaxy_table_stacking(fitspath,dataset, new_name)
    print('run complete')

def ind_detection(fitspath, dataset, bin_id):

    get_det3_tab = asc.read(fitspath+'get_det3_table.tbl')
    bin_tab = asc.read(fitspath+dataset+'_2d_binning_datadet3.tbl')
    N_gal_tab = asc.read(fitspath+dataset+'_Average_R23_O32_Values.tbl')
    stackmeas_tab = asc.read(fitspath+dataset+'_temperatures_metalicity.tbl')

    #From tables
    Source_id = get_det3_tab['Individual_IDs']
    O4959 = get_det3_tab['O4959']
    O5007 = get_det3_tab['O5007']
    Bin_number = bin_tab['Bin_number']
    O2 = get_det3_tab['O2']
    O3 = get_det3_tab['O3']
    Hb = get_det3_tab['Hb']
    N_Galaxies = N_gal_tab['N_stack']
    temp_bin = stackmeas_tab['T_e']
    
    R23 = get_det3_tab['R23']
    O32 = get_det3_tab['O32']
    
    #Initializing Arrays
    
    Source_IDs = []
    Bin_ID = []
    two_beta = []
    three_beta = []
    OIII4959 = []
    OIII5007=[]
    HBeta = []
    average_temp = []
    R23_ind = []
    O32_ind = []
    
    for ii in range(len(O2)):
        if Bin_number[ii] == bin_id:
            #print 'Bin_number:', Bin_number[ii], 'O2:', O2[ii], 'O3:', O3[ii], 'Hb:', Hb[ii]
            Source_IDs.append(Source_id[ii])
            Bin_ID.append(bin_id)
            R23_ind.append(R23[ii])
            O32_ind.append(O32[ii])
            two_beta.append(O2[ii]/Hb[ii])
            three_beta.append(O3[ii]/Hb[ii])
            OIII4959.append(O4959[ii])
            OIII5007.append(O5007[ii])
            HBeta.append(Hb[ii])
            average_temp.append(temp_bin[bin_id])
            
    
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(bin_id)+'_individual_ratios_temp.tbl'
    n = ('Source_ID','Bin_ID','Individual_R23', 'Individual_O32','two_beta', 'three_beta', 'OIII4959','OIII5007','HBeta','Temperature')   #'ID', 'R23_Average', 'O32_Average'
    ind_tab = Table([Source_IDs, Bin_ID, R23_ind, O32_ind, two_beta, three_beta, OIII4959, OIII5007, HBeta, average_temp], names=n) #ID, R23, O32,
    asc.write(ind_tab, individual_ascii, format = 'fixed_width_two_line')


def individual_galaxy_table_stacking(fitspath,dataset, new_name):
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    table_files = glob.glob(individual_ascii)
    table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        print asc_tab[0]
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,new_name, format='fixed_width_two_line', overwrite = True)





def individual_detection_MSC(fitspath, dataset, det3=True):
    '''
    Purpose: import all the required files to run composite_indv_detect.main from Metallicity Stack Commons 
    Out: ascii file: individual_derived_properties.tbl
    '''
    composite_file = fitspath +'bin_derived_properties.tbl'
    indv_em_line_file = fitspath + 'individual_properties.tbl'#file containing emission line information for individual galaxy

    indv_bin_file = fitspath + 'individual_bin_info.tbl' #bin information for each galaxy
    outfile = fitspath + 'individual_derived_properties.tbl'
    main(fitspath, dataset, composite_file, indv_em_line_file, indv_bin_file, outfile, det3=True)
    



