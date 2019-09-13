### Finds the temperature and metallicity for individual spectra in a given bin
##Version 2

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

from . import R_temp_calcul



def run_ind_detection(fitspath, dataset):    #, average_value_ascii):
    #N_gal_tab = asc.read(average_value_ascii)  #fitspath+dataset+'_Average_R23_O32_Values.tbl'
    ind_detection(fitspath,dataset)
    new_name = fitspath+'Individual_ratio_temperature.tbl'
    individual_galaxy_table_stacking(fitspath,dataset, new_name)
    print('run complete')

def ind_detection(fitspath, dataset):

    get_det3_tab = asc.read(fitspath+'get_det3_table2.tbl')
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
    N_Galaxies = N_gal_tab['N_Galaxies']
    Bin_ID = N_gal_tab['ID']
    temp_bin = stackmeas_tab['Temperature']
    
    R23 = get_det3_tab['R23']
    O32 = get_det3_tab['O32']
      
    for ii in Bin_ID:
        current_bin = np.where((Bin_number == ii))[0]
        Source_IDs = Source_id[current_bin]
        Bin_ID = np.repeat(ii, len(current_bin))
        two_beta = (O2[current_bin]/Hb[current_bin])
        three_beta = (O3[current_bin]/Hb[current_bin])
        OIII4959 = O4959[current_bin]
        OIII5007=O5007[current_bin]
        HBeta = Hb[current_bin]
        average_temp = np.repeat(temp_bin[ii],len(current_bin))  #np.repeat() or [ii]* number of times to repeat
        R23_ind = R23[current_bin]
        O32_ind = O32[current_bin]


        individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(ii)+'_individual_ratios_temp.tbl'
        n = ('Source_ID','Bin_ID','Individual_R23', 'Individual_O32','two_beta', 'three_beta', 'OIII4959','OIII5007','HBeta','Temperature')   #'ID', 'R23_Average', 'O32_Average'
        ind_tab = Table([Source_IDs, Bin_ID, R23_ind, O32_ind, two_beta, three_beta, OIII4959, OIII5007, HBeta, average_temp], names=n) #ID, R23, O32,
        asc.write(ind_tab, individual_ascii, format = 'fixed_width_two_line')


def individual_galaxy_table_stacking(fitspath,dataset, new_name):
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    table_files = glob.glob(individual_ascii)
    #table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        print asc_tab[0]
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,new_name, format='fixed_width_two_line', overwrite = True)






def call_metallicity_calculation(fitspath, dataset, individual_combine_table):
    #individual_combine_table = /Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_ratio_temperature.tbl'

    
    out_ascii = fitspath+'Individual_temp_metal.tbl'
    out_fits = fitspath+'Individual_temp_metal.fits'
    pdf_name = fitspath+'Individual_temp_metal.pdf'
    dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'

    R_temp_calcul.run_function(fitspath, dataset, out_ascii, out_fits, pdf_name, individual_combine_table, dust_ascii, dustatt=False)

    
    Metal_Te_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_temp_metal.tbl'
    MT_ascii = asc.read(Metal_Te_table)

    R23_ind = MT_ascii['R23']
    O32_ind = MT_ascii['O32']
    ind_metal = MT_ascii['com_O_log']
    R23 = np.log10(R23_ind)
    O32 = np.log10(O32_ind)

    pdf_name= 'R23vsO32_individual_metal_plots.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    plot1= ax1.scatter(R23,O32,0.8, c=ind_metal ,marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('R23')
    ax1.set_ylabel('O32')
    ax1.set_title('R23 vs. O32  Colormap=Metallicity')
    fig1.set_size_inches(8,8)
    fig1.savefig(pdf_pages, format='pdf')

    pdf_pages.close()


def massvslum_plots(fitspath, dataset, individual_combine_table):
    out_ascii = fitspath+'Individual_temp_metal.tbl'
    out_fits = fitspath+'Individual_temp_metal.fits'
    pdf_name = fitspath+'Individual_temp_metal.pdf'
    dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'

    R_temp_calcul.run_function(fitspath, dataset, out_ascii, out_fits, pdf_name, individual_combine_table, dust_ascii, dustatt=False)

    
    Metal_Te_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_temp_metal.tbl'
    MT_ascii = asc.read(Metal_Te_table)

    MASS = MT_ascii['Mass']
    LHBETA = MT_ascii['LHbeta']
    ind_metal = MT_ascii['com_O_log']

    pdf_name= 'MvsL_individual_metal_plots.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    plot1= ax1.scatter(MASS,LHBETA,0.8, c=ind_metal ,marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('Mass')
    ax1.set_ylabel('LHBETA')
    ax1.set_title('Mass vs. Luminosity Colormap=Metallicity')
    fig1.set_size_inches(8,8)
    fig1.savefig(pdf_pages, format='pdf')
    pdf_pages.close()
