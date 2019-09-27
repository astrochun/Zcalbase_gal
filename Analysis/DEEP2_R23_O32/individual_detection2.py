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

#Caroline's where statement:  valid_ind = np.where((np.isfinite(temp_x)==True) & (temp_x>0) & (temp_x < 1e13) & (temp_flag == 0))[0]

def individual_galaxy_detections_plots(fitspath,dataset, new_name):
    #fitspath+dataset+'_Average_R23_O32_Values.tbl'
    #new_name = fitspath+'Individual_ratio_temperature.tbl'

    #Rewrite this function so that it pulls values for both Caroline and Reagen's individual detections
    ind_detection(fitspath,dataset)
    
    ###Write code so that you pulls the new table, don't think we need to stack them anymore
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    '''table_files = glob.glob(individual_ascii)
    table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        print asc_tab[0]
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,new_name, format='fixed_width_two_line', overwrite = True)'''


    pdf_pages= PdfPages()
    if: ####Something that calls Reagen's data:
        individual_combine_table = ''
        #out_ascii = fitspath+'Individual_temp_metal.tbl'
        #out_fits = fitspath+'Individual_temp_metal.fits'
        #temp_pdf = fitspath+'Individual_temp_metal.pdf'
        #dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'
        metal_pdf= fitspath + 'R23vsO32_individual_metal_plots.pdf'
        call_metallicity_calculation(fitspath,dataset,individual_combine_table,out_ascii,out_fits,temp_pdf, metal_pdf)  #verification table?!


    if: ####Something that calls Caroline's data:
        individual_combine_table = ''
        #out_ascii = fitspath+'Individual_temp_metal.tbl'
        #out_fits = fitspath+'Individual_temp_metal.fits'
        #temp_pdf = fitspath+'Individual_temp_metal.pdf'
        #dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'
        metal_pdf= fitspath + 'MassvsLum_individual_metal_plots.pdf'
        massvslum_plots(fitspath, dataset, individual_combine_table,out_ascii,out_fits, temp_pdf, metal_pdf)  #verification table?!

    
def ind_detection(fitspath, dataset):

    get_det3_tab = asc.read(fitspath+'get_det3_table2.tbl')
    get_massL_tab = asc.read()
    
    


    bin_tab = asc.read(fitspath+dataset+'_2d_binning_datadet3.tbl')
    N_gal_tab = asc.read(fitspath+dataset+'_Average_R23_O32_Values.tbl')
    stackmeas_tab = asc.read(fitspath+dataset+'_temperatures_metalicity.tbl')

    #From tables
    Source_id = get_det3_tab['Individual_IDs']
    O4959 = get_det3_tab['O4959']
    O5007 = get_det3_tab['O5007']
    Hgamma = get_det3_tab['Hgamma']
    SNRHG = get_det3_tab['SNRHG']
    SNR4363 = get_det3_tab['SNR4363']
    O4363 = get_det3_tab['O4363']
    Bin_number = bin_tab['Bin_number']
    O2 = get_det3_tab['O2']
    O3 = get_det3_tab['O3']
    Hb = get_det3_tab['Hb']
    N_Galaxies = N_gal_tab['N_Galaxies']
    Bin_ID = N_gal_tab['ID']
    temp_bin = stackmeas_tab['Temperature']
    detection = stackmeas_tab['Detection']
    
    R23 = get_det3_tab['R23']
    O32 = get_det3_tab['O32']
    
    for ii in Bin_ID:
        current_bin = np.where((Bin_number == ii))[0]
        Source_IDs = Source_id[current_bin]
        Bin_ID = np.repeat(ii, len(current_bin))
        two_beta = (O2[current_bin]/Hb[current_bin])
        three_beta = (O3[current_bin]/Hb[current_bin])
        OIII4959 = O4959[current_bin]
        OIII5007 = O5007[current_bin]
        HBeta = Hb[current_bin]
        average_temp = np.repeat(temp_bin[ii],len(current_bin))  #np.repeat() or [ii]* number of times to repeat
        Detection = np.repeat(detection[ii],len(current_bin))
        R23_ind = R23[current_bin]
        O32_ind = O32[current_bin]

        OIII4363 = O4363[current_bin]
        HGamma = Hgamma[current_bin]
        SNR_HG = SNRHG[current_bin]
        SNR_4363 = SNR4363[current_bin]

        Reagens_data = np.where
        Carolines_data = np.where
        


        individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(ii)+'_individual_ratios_temp.tbl'
        n = ('Source_ID','Bin_ID','Individual_R23', 'Individual_O32','two_beta', 'three_beta', 'HGamma','SNR_HG','OIII4363','SNR_4363','OIII4959','OIII5007','HBeta','Temperature','Detection')   #'ID', 'R23_Average', 'O32_Average'
        ind_tab = Table([Source_IDs, Bin_ID, R23_ind, O32_ind, two_beta, three_beta, HGamma, SNR_HG, OIII4363, SNR_4363, OIII4959, OIII5007, HBeta, average_temp,Detection], names=n) #ID, R23, O32,
        asc.write(ind_tab, individual_ascii, format = 'fixed_width_two_line')







def call_metallicity_calculation(fitspath, dataset, individual_combine_table):
    #individual_combine_table = /Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_ratio_temperature.tbl'
    #out_ascii = fitspath+'Individual_temp_metal.tbl'
    #out_fits = fitspath+'Individual_temp_metal.fits'
    #pdf_name = fitspath+'Individual_temp_metal.pdf'
    #dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'

    R_temp_calcul.run_function(fitspath, dataset, out_ascii, out_fits, pdf_name, individual_combine_table, dust_ascii, dustatt=False)

    
    Metal_Te_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_temp_metal.tbl'
    MT_ascii = asc.read(Metal_Te_table)

    R23_ind = MT_ascii['R23']
    O32_ind = MT_ascii['O32']
    Detection = MT_ascii['Detection']
    ind_metal = MT_ascii['com_O_log']
    R23_log = np.log10(R23_ind)
    O32_log = np.log10(O32_ind)
    detect = np.where((Detection ==1))[0]
    nan_detect = np.where((Detection ==0))[0]
    
    pdf_pages = PdfPages(metal_pdf)

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    plot1= ax1.scatter(R23[detect],O32[detect],0.8, c=ind_metal ,marker='*')
    plot2= ax1.scatter(R23[nan_detect],O32[nan_detect],0.8, c=ind_metal ,marker='.')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('R23')
    ax1.set_ylabel('O32')
    ax1.set_title('R23 vs. O32  Colormap=Metallicity')
    fig1.set_size_inches(8,8)
    fig1.savefig(pdf_pages, format='pdf')

    pdf_pages.close()


def massvslum_plots(fitspath, dataset, individual_combine_table,out_ascii,out_fits, temp_pdf, metal_pdf): 
    #out_ascii = fitspath+'Individual_temp_metal.tbl'
    #out_fits = fitspath+'Individual_temp_metal.fits'
    #temp_pdf = fitspath+'Individual_temp_metal.pdf'
    #dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'

    R_temp_calcul.run_function(fitspath, dataset, out_ascii, out_fits, pdf_name, individual_combine_table, dust_ascii, dustatt=False)

    
    Metal_Te_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_temp_metal.tbl'
    MT_ascii = asc.read(Metal_Te_table)

    MASS = MT_ascii['Mass']
    LHBETA = MT_ascii['LHbeta']
    ind_metal = MT_ascii['com_O_log']

    detect = np.where((LHBETA =! -1.0))[0]

    pdf_name= 'MvsL_individual_metal_plots.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    plot1= ax1.scatter(MASS,LHBETA[detect],0.8, c=ind_metal ,marker='*')
    cb = fig1.colorbar(plot1)
    cb.set_label('Metallicity')
    ax1.set_xlabel('Mass')
    ax1.set_ylabel('LHBETA')
    ax1.set_title('Mass vs. Luminosity Colormap=Metallicity')
    fig1.set_size_inches(8,8)
    fig1.savefig(pdf_pages, format='pdf')
    pdf_pages.close()
