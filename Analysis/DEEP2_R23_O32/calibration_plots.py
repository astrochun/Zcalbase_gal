###Calibrating my data to other work 
###Question: Should we organize into the nan dections and detections before making these plots???

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
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ
import glob

from Zcalbase_gal.Analysis import local_analog_calibration, green_peas_calibration

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/'
dataset  = 'n_Bins'
temp_tab = fitspath+ 'bin_derived_properties.tbl'
verification_table = fitspath + 'bin_validation.revised.tbl'  

def LAC_GPC_plots():           #fitspath, dataset,temp_tab,verification_table):
    
    out_pdf = fitspath+ '/'+dataset+'_LAC.pdf'

    valid = asc.read(verification_table)
    detect = valid['Detection']
    print(detect)
    
    temp_table= asc.read(temp_tab)
    SN_4363 = temp_table['OIII_4363_S/N']
    det_4363 = np.where((detect == 1.0))[0]
    print('det_4363: ', det_4363)
    rlimit = np.where((detect == 0.5))[0]
    print('Begin Local analog Calibration')
    
    ###Implimenting Local analog calibration###
    
    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')
    
    #DEEP2 Derived 
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_OH = derived['OH'].data
    
    #MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_MACT['OH'].data
    
    
    
    O32_all = temp_table['logO32_Composite']
    print(O32_all)
    R23_all = temp_table['logR23_Composite']
    com_O_log = temp_table['12+log(O/H)']  #This is the 12+log(OH) value
    ID = temp_table['bin_ID']
    #print O32_all
    
    det_O32 = O32_all[det_4363]
    print('det_O32: ', det_O32)
    det_R23 = R23_all[det_4363]
    det_OH  = com_O_log[det_4363]
    det_ID  = ID[det_4363]
    
    nandet_O32 = O32_all[rlimit]
    nandet_R23 = R23_all[rlimit]
    #print(len(nandet_O32), len(nandet_R23))
    nandet_OH  = com_O_log[rlimit]
    nandet_ID  = ID[rlimit]

    #Individual Detections from Zcalbase_gal Analysis
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/individual_derived_properties.tbl'     
    individual = asc.read(individual_ascii)
    logR23 = individual['logR23']
    logO32 = individual['logO32']
    com_log = individual['12+log(O/H)']
    bin_ID = individual['bin_ID']


    pea_out_pdf1 = fitspath+ '/'+dataset+'_AGPC_valid.pdf'
    pea_out_pdf2 = fitspath+ '/'+dataset+'_GPC_limits.pdf'
    pea_out_pdf3 = fitspath+ '/'+dataset+'_GPC_zcalbase_all.pdf'
    out_pdf_LAC = fitspath+ '/'+ dataset+'ALAC_plot.pdf'
    
    label = ['Detection','Robust Limits','DEEP2', 'MACT']
    marker = ['D',r'$\uparrow$','3','4']
    
    rlR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
    rlO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
    rOH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]

    green_peas_calibration.main(rlR23,rlO32, rOH, pea_out_pdf2, n_bins=6, xra=[0.3,1.15], yra=[6.5,9.10], marker=marker, edgecolors= ['face','face', 'none', 'none'], alpha = [0.5, 0.5, 0.5, 0.5], label=label, fit=False, silent=False, verbose=True)
    # marker=['.','*','^','o'], label=['Detection','Non-Dectection','DEEP2', 'MACT']

    
    lR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
    lO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
    OH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]
    
    local_analog_calibration.main(lR23, lO32, OH, out_pdf_LAC, yra=[7.0,9.0], ctype=['b','g','r','m'], label=label, marker= marker, silent=False, verbose=True)
    print('finished LAC plot') 


    lR23 = [det_R23,der_R23,der_R23_MACT]
    lO32 = [det_O32,der_O32,der_O32_MACT]
    OH   = [det_OH, der_OH, der_OH_MACT]
    
    #green_peas_calibration.main(lR23,lO32, OH, pea_out_pdf1, n_bins=6, xra=[0.3,1.15], yra=[6.5,9.10], marker=['D','3','4'], edgecolors= ['face','face', 'none'], alpha = [0.5, 0.5, 0.5], label=['Detection','DEEP2', 'MACT'], fit=False, silent=False, verbose=True)








'''
    lR23 = [det_R23,der_R23,der_R23_MACT]
    
    lO32 = [det_O32,der_O32,der_O32_MACT]
    
    OH   = [det_OH, der_OH, der_OH_MACT]
    
    green_peas_calibration.main(lR23,lO32, OH, pea_out_pdf1, n_bins=6, xra=[0.3,1.15], yra=[6.5,9.10], marker=['D','3','4'], edgecolors= ['face','face', 'none'], label=['Detection','DEEP2', 'MACT'], fit=False, silent=False, verbose=True)
    # marker=['.','^','o'], label=['Detection','Non-Dectection','DEEP2', 'MACT']
    print('Done with detections.')

    
    rlR23 = [det_R23,nandet_R23, logR23]
    rlO32 = [det_O32,nandet_O32, logO32]
    rOH   = [det_OH,nandet_OH, com_log] 

    #green_peas_calibration.main(rlR23,rlO32, rOH, pea_out_pdf3, n_bins=6, xra=[0.3,1.15], yra=[6.5,9.10], marker=['D','X', '.'],edgecolors= ['face','face', 'none'], alpha = [1, 1, 0.2], label=['Detection','Robust Limits', 'Zcalbase_gal'], fit=False, silent=False, verbose=True)
    



        
    if dataset == 'R23_Grid':
        lR23 = [det_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,der_O32,der_O32_MACT]
        OH   = [det_OH, der_OH, der_OH_MACT]
        local_analog_calibration.main(lR23, lO32, OH, out_pdf, ctype=['b','r','m'], label=['Detection','DEEP2', 'MACT'], silent=False, verbose=True)
        
    if dataset == 'O32_Grid' or dataset == 'Grid':    
        lR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
        OH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]
        local_analog_calibration.main(lR23, lO32, OH, out_pdf, ctype=['b','g','r','m'], label=['Detection','Non-Dectection','DEEP2', 'MACT'], silent=False, verbose=True)
    # ID=[det_ID,nandet_ID]

    if dataset == 'Voronoi10' or dataset == 'Voronoi14' or dataset == 'Voronoi20' or dataset =='Double_Bin':
        lR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
        OH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]

        local_analog_calibration.main(lR23, lO32, OH, out_pdf, yra=[7.0,9.0], ctype=['b','g','r','m'], label=['Detection','Non-Dectection','DEEP2','MACT'], silent=False, verbose=True)

    if dataset == 'n_Bins':
        lR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
        OH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]

        #local_analog_calibration.main(lR23, lO32, OH, out_pdf, yra=[7.0,9.0], ctype=['b','g','r','m'], label=['Detection','Non-Dectection','DEEP2','MACT'], silent=False, verbose=True)
        print('finished LAC plot') 








    
    ###Green Pea Calibration###
    pea_out_pdf = fitspath+ '/'+dataset+'_GPC.pdf'

    if dataset == 'R23_Grid':
        lR23 = [det_R23,der_R23,der_R23_MACT]                    #[det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,der_O32,der_O32_MACT]                    #[det_O32,nandet_O32,der_O32,der_O32_MACT]
        OH   = [det_OH, der_OH, der_OH_MACT]                    #[det_OH,nandet_OH, der_OH, der_OH_MACT]

    else: 
        lR23 = [det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32,nandet_O32,der_O32,der_O32_MACT]
        print('lO32:', lO32)
        OH   = [det_OH,nandet_OH, der_OH, der_OH_MACT]'''



def individual_GPC(individual_ascii, validation_table):

    pea_out_pdf_ind = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/jiang_plot_individual.pdf'
    individual = asc.read(individual_ascii)
    logR23 = individual['logR23']
    logO32 = individual['logO32']
    com_log = individual['12+log(O/H)']
    bin_ID = individual['bin_ID']

    valid = asc.read(validation_table)
    Detections = valid['Detection']
    detect = np.where((Detections == 1.0))[0]
    rlimit = np.where((Detections == 0.5))[0]
    bins = valid['bin_ID']

    ID_detect = bins[detect]
    ID_rlimit = bins[rlimit]


    lR23 = [logR23]
    lO32 = [logO32]
    OH   = [com_log]

    green_peas_calibration.main(lR23,lO32, OH, pea_out_pdf_ind, n_bins=6, xra=[0.3,1.15], yra=[6.5,9.10], marker=['3'], label=['Individual Detection'], fit=False, silent=False, verbose=True)

    
        
   
