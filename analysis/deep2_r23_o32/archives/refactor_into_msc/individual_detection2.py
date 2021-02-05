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
from . import general, get_det3

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

def individual_galaxy_detections_plots(fitspath,project, dataset):
    #fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/practice/'

    ###Get the indexing statements
    get_det3(fitspath, individual_detect = True)
    

    ###Import Table
    #R23O32_tab = asc.read(fitspath+'get_det3_table.tbl')
    #ML_tab = asc.read(fitspath+'mass_LHbeat_tab.tbl')
    #where_tab = asc.read(fitspath+'indexed_individual.tbl')


    #Rewrite this function so that it pulls values for both Caroline and Reagen's individual detections
    ind_detection(fitspath,project,dataset)
    
    


    #pdf_pages= PdfPages()
    if project == 'R23O32': 
        new_name = 'R23O32_individual_detect.tbl'
        stacking_ascii(fitspath, dataset, new_name)
        individual_combine_table = fitspath+new_name
        out_ascii = fitspath+'Individual_temp_metal.tbl'
        out_fits = fitspath+'Individual_temp_metal.fits'
        temp_pdf = fitspath+'Individual_temp_metal.pdf'
        dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'
        metal_pdf= fitspath + 'R23vsO32_individual_metal_plots.pdf'
        call_metallicity_calculation(fitspath,dataset,individual_combine_table,out_ascii,out_fits,temp_pdf, metal_pdf, dust_ascii)  #verification table?!


    if project == 'ML':
        individual_combine_table = ''
        out_ascii = fitspath+'Individual_temp_metal.tbl'
        out_fits = fitspath+'Individual_temp_metal.fits'
        temp_pdf = fitspath+'Individual_temp_metal.pdf'
        dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'
        metal_pdf= fitspath + 'MassvsLum_individual_metal_plots.pdf'
        massvslum_plots(fitspath, dataset, individual_combine_table,out_ascii,out_fits, temp_pdf, metal_pdf, dust_ascii)  #verification table?!

    
def ind_detection(fitspath,project,dataset):
    get_det3_tab= asc.read(fitspath+'get_det3_table.tbl')
    ml_tab = asc.read(fitspath+'mass_LHbeta_tab.tbl')
    
    if project =='R23O32':
        bin_tab = asc.read(fitspath+dataset+'_2d_binning_datadet3.tbl')
        N_gal_tab = asc.read(fitspath+dataset+'_Average_R23_O32_Values.tbl')
        stackmeas_tab = asc.read(fitspath+dataset+'_temperatures_metalicity.tbl')
        
        #From tables
        Source_id = get_det3_tab['Individual_IDs']
        O4959 = get_det3_tab['O4959']
        O5007 = get_det3_tab['O5007']
        Hgamma = get_det3_tab['HGAMMA']
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
        
            individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/'+str(ii)+'_individual_ratios_temp.tbl'
        
            n = ('Source_ID','Bin_ID','Individual_R23', 'Individual_O32','two_beta', 'three_beta', 'HGamma','SNR_HG','OIII4363','SNR_4363','OIII4959','OIII5007','HBeta','Temperature','Detection')   #'ID', 'R23_Average', 'O32_Average'
            ind_tab = Table([Source_IDs, Bin_ID, R23_ind, O32_ind, two_beta, three_beta, HGamma, SNR_HG, OIII4363, SNR_4363, OIII4959, OIII5007, HBeta, average_temp,Detection], names=n) #ID, R23, O32,
            asc.write(ind_tab, individual_ascii, format = 'fixed_width_two_line')















    if project == 'ML_data':
        ###INPUT FILES 
        #contains bin indices
        HB_bin_file = fitspath_ini+'From_Caroline/mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz'
        
        #contains electron temperatures
        HB_derived_prop_file = fitspath_ini+'From_Caroline/hbeta_revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
        mass_derived_prop_file = fitspath_ini+'From_Caroline/revised_75_112_113_300_600_1444_1444_updated_massbin_derived_properties_metallicity.tbl'
        
        #contains mass and luminosity values
        mass_LHbeta_file = fitspath_ini+'From_Caroline/revised_75_112_113_300_600_1444_1444_mass_SFR_data.npz'
        
        #contains combined visual detections and S/N > 3 detections
        HB_valid_file = fitspath_ini+'From_Caroline/hbeta_revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'
        mass_valid_file = fitspath_ini+'From_Caroline/revised_75_112_113_300_600_1444_1444_massbin_validation.tbl'


        HB_bin_info = np.load(path2 + HB_bin_file)
        HB_bin_Te = asc.read(path2 + HB_derived_prop_file)
        mass_bin_Te = asc.read(path2 + mass_derived_prop_file)
        mass_LHbeta = np.load(path2 + mass_LHbeta_file)
        HB_valid_table = asc.read(path2 + HB_valid_file)
        mass_valid_table = asc.read(path2 + mass_valid_file)
        
        
        HB_bin_ind = HB_bin_info['bin_ind']                   #valid HBeta bin indices relative to unsorted data table
        HB_Te = HB_bin_Te['Temperature'].data                 #HBeta bin electron temperatures
        mass_Te = mass_bin_Te['Temperature'].data             #Mass bin electron temperatures
        mass = mass_LHbeta['mass']
        LHbeta = mass_LHbeta['lum']
        HB_detections = HB_valid_table['Detection'].data      #HBeta bin detections
        mass_detections = mass_valid_table['Detection'].data  #Mass bin detections
        
        hdu = fits.open(path3 + 'DEEP2_all_line_fit.fits')
        line_table = hdu[1].data                  
        
        
        #HBeta bin valid indices
        HB_valid_idx = []
        for ii in range(len(HB_bin_ind)):                      #creates 1D array of all indices
            HB_valid_idx += list(HB_bin_ind[ii])               
        idx_array = np.array(HB_valid_idx)                     #len = 4088
            
        EBV_array = np.zeros(len(idx_array))             
        indiv_detect_array = np.zeros(len(idx_array))        
        mass_array = np.log10(mass[idx_array])
        LHbeta_array = LHbeta[idx_array]        
    
        line_table = line_table[idx_array]              
        objno = line_table['OBJNO']
        OII_Flux = line_table['OII_FLUX_DATA']
        OII_SN = line_table['OII_SNR']
        OIII4959_Flux = line_table['OIIIB_FLUX_DATA']
        OIII4959_SN = line_table['OIIIB_SNR']
        OIII5007_Flux = line_table['OIIIR_FLUX_DATA']
        OIII5007_SN = line_table['OIIIR_SNR']
        HBETA_Flux = line_table['HB_FLUX_DATA']
        HBETA_SN = line_table['HB_SNR'] 
    
        #only works if mass bins are split into two SFR bins --> fix library.py and adapt for more general case
        mass_bin_ID = []
        LHbeta_bin_ID = []
        HB_Te_array = []
        mass_Te_array = []
        HB_detect_array = []
        mass_detect_array = []
        kk = 1    
        for ii in range(len(HB_bin_ind)):
            for jj in range(len(HB_bin_ind[ii])):
                mass_bin_ID.append(kk)
                LHbeta_bin_ID.append(ii + 1)
                HB_Te_array.append(HB_Te[ii])
                mass_Te_array.append(mass_Te[kk - 1])
        
                #Excluding the last bin --> S/N > 3, but possible AGN contribution
                if HB_bin_type == 'hbeta_revised_75_112_113_300_600_1444_1444':
                    if (ii + 1 == 14):
                        HB_detect_array.append(0.0)
                    else:
                        HB_detect_array.append(HB_detections[ii])    
                
                if mass_bin_type == 'revised_75_112_113_300_600_1444_1444':
                    mass_detect_array.append(mass_detections[kk - 1])
            
            if (ii + 1) % 2 == 0:
                kk += 1
    
        out_ascii = path2 + 'indivgals_Te_lineRatios.tbl'
        n = ('OBJNO', 'Mass_Bin_ID', 'LHBeta_Bin_ID', 'Log10(Mass)', 'HBeta_Luminosity', 'Mass Bin Te', 'LHBeta Bin Te', 'OII_Flux', 'OII_SN', 'OIII4959_Flux', 'OIII4959_SN', 'OIII5007_Flux', 'OIII5007_SN', 'HBETA_Flux', 'HBETA_SN', 'Mass Bin Detections', 'LHBeta Bin Detections', 'Individual Detections', 'EBV')
        Te_ratio_table = Table([objno, mass_bin_ID, LHbeta_bin_ID, mass_array, LHbeta_array, mass_Te_array, HB_Te_array,OII_Flux, OII_SN, OIII4959_Flux, OIII4959_SN, OIII5007_Flux, OIII5007_SN, HBETA_Flux,HBETA_SN, mass_detect_array, HB_detect_array, indiv_detect_array, EBV_array], names = n)

        asc.write(Te_ratio_table, out_ascii, format = 'fixed_width_two_line', overwrite = True)
    
    
        hdu.close()

    if project == '': 
        all_individual_values = fitspath+'indexed_individual.tbl'
        all_ind = asc.read('all_individual_values')



def stacking_ascii(fitspath, dataset,new_name):
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    table_files = glob.glob(individual_ascii)
    #table_files.sort()
    print('Stacking tables')
    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        #print asc_tab[0]
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,fitspath+new_name, format='fixed_width_two_line', overwrite = True)



def call_metallicity_calculation(fitspath, dataset, individual_combine_table, out_ascii,out_fits, temp_pdf, metal_pdf, dust_ascii):
    #individual_combine_table = /Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_ratio_temperature.tbl'
    #out_ascii = fitspath+'Individual_temp_metal.tbl'
    #out_fits = fitspath+'Individual_temp_metal.fits'
    #pdf_name = fitspath+'Individual_temp_metal.pdf'
    #dust_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'

    R_temp_calcul.run_function(fitspath, dataset, out_ascii, out_fits, temp_pdf, individual_combine_table, dust_ascii, dustatt=False)

    
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

    detect = np.where((LHBETA != -1.0))[0]

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








'''Not In Use

###Write code so that you pull the new table, don't think we need to stack them anymore
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    table_files = glob.glob(individual_ascii)
    table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        print asc_tab[0]
        if ii == 0: vstacking = asc_tab
        else: vstacking = vstack([vstacking,asc_tab])
    asc.write(vstacking,new_name, format='fixed_width_two_line', overwrite = True)





'''
