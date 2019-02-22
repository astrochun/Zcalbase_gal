#Writes Ascii Files with all the Correct data in them

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


###Verified Sources###
ver_20 = np.array([])
ver_14 = np.array([0,7,10,11,12])
ver_10 = np.array([1,9,18,21])
ver_Grid = np.array([11,13,19,20,21])
ver_sR23 = np.array([0,4,5,6])
ver_sO32 = np.array([6])





'''fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi20_0127/'
dataset = 'Voronoi20'
temp_table = fitspath+ '/'+dataset+'_temperatures_metalicity.tbl'
combine_flux_ascii = fitspath+dataset+'_combined_flux_table.tbl'
temp_tab = asc.read(temp_table)
com_flux = asc.read(combine_flux_ascii)
ver_tab = fitspath+'/'+dataset+'_verification.tbl'''


def verification(fitspath, dataset,temp_table, combine_flux_ascii, ver_tab): 
    temp_tab = asc.read(temp_table)
    com_flux = asc.read(combine_flux_ascii)

    SN_4363 = temp_tab['S/N_4363'].data
    O32_all = temp_tab['O32_Composite'].data
    R23_all = temp_tab['R23_Composite'].data
    com_O_log = temp_tab['com_O_log']  #This is the 12+log(OH) value
    ID = temp_tab['ID'].data
    all_4363 =temp_tab['Observed_Flux_4363'].data
    N_Gal_all = temp_tab['N_Galaxies'].data
    temp_all = temp_tab['Temperature'].data

    #Detection Where Statements
    if dataset== 'Voronoi20': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'Voronoi14': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'Voronoi10': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'Grid': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'R23_Grid': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'O32_Grid': det_4363 = np.where((SN_4363>=3))[0]
    if dataset== 'Double_Bin': det_4363 = np.where((ID ==0)| (ID ==1)| (ID ==2) | (ID ==7) | (ID ==9) | (ID ==10) | (ID ==11) | (ID ==13))[0]
    #if dataset== 'Double_Bin': det_4363 = np.where((SN_4363>=3) | (ID ==0)| (ID ==1)| (ID ==2) | (ID ==7) | (ID ==9) | (ID ==10) | (ID ==11) | (ID ==13))[0]
    print 'det_4363:', det_4363
    
   

    det_O32 = O32_all[det_4363]
    det_R23 = R23_all[det_4363]
    det_OH  = com_O_log[det_4363]
    det_ID  = ID[det_4363]
    det_SN_4363 = SN_4363[det_4363]
    det_obs4363 = all_4363[det_4363]
    det_N_Gal = N_Gal_all[det_4363]
    det_temp = temp_all[det_4363]
    
    
    n =('ID', '4363', 'S/N','R23', 'O32', 'N_Galaxies','Temperature','Com_O_log')
    tab1 = Table([det_ID, det_obs4363,det_SN_4363,det_R23, det_O32,det_N_Gal,det_temp,det_OH], names=n)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')



###Equivalent Width correction
'''bottom side of the iceburg / continum 
take negative component of gauss and subtract off the postive component 

total gauss - the positve component? 
x = double gaussian of o1[0,1,2] = 0 
x-o1[3] from [-2.5sigma to 2.5 sigma]
equivalent width in terms of angstrumns 
update plots '''
