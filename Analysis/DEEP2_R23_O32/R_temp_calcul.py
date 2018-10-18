#Calculates the R value, electron temperature, and metallicity from the flux table
#produced by the zoom_and_gauss_general functions
#
#Currently running: Grid
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
from mpl_toolkits.mplot3d import Axes3D
import sys

#For generalizing for several users
from getpass import getuser

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
'''
#Voronoi10
#Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi10/Voronoi_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Grid
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/Grid_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Voronoi14
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Ascii Table Calls 
OIII5007 = data1['OIII_5007_Flux_Observed'].data
OIII4959 = data1['OIII_4958_Flux_Observed'].data
OIII4363 = data1['OIII_4363_Flux_Observed'].data
HBETA    = data1['HBETA_Flux_Observed'].data
OII3727  = data1['OII_3727_Flux_Observed'].data
R23_avg      = data1['R_23_Average'].data
O32_avg      = data1['O_32_Average'].data
N_Galaxy = data1['N_Galaxies'].data

SN_5007       = data1['OIII_5007_S/N'].data
SN_4959       = data1['OIII_4958_S/N'].data
SN_4363       = data1['OIII_4363_S/N'].data
SN_HBETA      = data1['HBETA_S/N'].data
SN_3727       = data1['OII_3727_S/N'].data


#Fits Table Calls (This is incorrect)
OIII5007 = header['OIII_5007_Flux_Observed']
OIII4959 =  header['OIII_4958_Flux_Observed']
OIII4363 =  header['OIII_4363_Flux_Observed']
HBETA    =  header['HBETA_Flux_Observed']
OII3727  =  header['OII_3727_Flux_Observed']
R23_avg      =  header['R_23_Average']
O32_avg      =  header['O_32_Average']
N_Galaxy =  header['N_Galaxies']

SN_5007       =  header['OIII_5007_S/N']
SN_4959       =  header['OIII_4958_S/N']
SN_4363       =  header['OIII_4363_S/N']
SN_HBETA      =  header['HBETA_S/N']
SN_3727       =  header['OII_3727_S/N']

R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
O32_composite = np.log10((1.33*OIII5007)/OII3727)
'''

#Constants

a = 13205
b = 0.92506
c = 0.98062
#x = S/N([O iii]4363)

#Calculations on page 37 of paper

def R_calculation(OIII4363, OIII5007, OIII4959):
  
    R_value = OIII4363/(OIII5007+OIII4959)
    return R_value

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)      #np.zeros(len(OIII5007))
  
    return T_e


def metalicity_calculation(T_e,OIII5007, OIII4959, OIII4363, HBETA, OII3727):
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17
    #What is x?
    '''t = 1e-4 * T_e
    x = 1e-4 * 1e3 * t**(-0.5)'''
    
    two_beta = OII3727/HBETA
    three_beta= (OIII4959+OIII5007)/HBETA
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


def run_function(fitspath, out_ascii, out_fits, pdf_name,  combine_flux_ascii):  #combine_fits, header
    #Fits Table Calls
    #combine_fits, header = fits.getdata(combine_flux_table, header = True)
    combine_fits= asc.read(combine_flux_ascii)

    #Ascii Table Calls 
    OIII5007 = combine_fits['OIII_5007_Flux_Observed'].data
    OIII4959 = combine_fits['OIII_4958_Flux_Observed'].data
    OIII4363 = combine_fits['OIII_4363_Flux_Observed'].data
    HBETA    = combine_fits['HBETA_Flux_Observed'].data
    OII3727  = combine_fits['OII_3727_Flux_Observed'].data
    R23_avg      = combine_fits['R_23_Average'].data
    O32_avg      = combine_fits['O_32_Average'].data
    N_Galaxy = combine_fits['N_Galaxies'].data

    SN_5007       = combine_fits['OIII_5007_S/N'].data
    SN_4959       = combine_fits['OIII_4958_S/N'].data
    SN_4363       = combine_fits['OIII_4363_S/N'].data
    SN_HBETA      = combine_fits['HBETA_S/N'].data
    SN_3727       = combine_fits['OII_3727_S/N'].data

    R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
    O32_composite = np.log10((1.33*OIII5007)/OII3727)


    R_value = R_calculation(OIII4363, OIII5007, OIII4959)
    T_e = temp_calculation(R_value)
    O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metalicity_calculation(T_e,OIII5007, OIII4959, OIII4363, HBETA, OII3727)
    
    print R_value
    print T_e
    print min(T_e), max(T_e)
    print O_s_ion, O_d_ion

   

    #Ascii Table
    #out_ascii = fitspath+ '/Grid_temperatures_metalicity_asc_table.tbl'
    #out_fits = fitspath+ '/Grid_temperatures_metalicity_asc_table.fits'
    if not exists(out_ascii):
        n=  ('R23_Composite', 'O32_Composite', 'R_23_Average', 'O_32_Average', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959', 'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727', 'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log')
        
        tab0 = Table([ R23_composite, O32_composite, R23_avg, O32_avg, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, T_e, O_s_ion, O_d_ion, com_O_log], names=n)
        asc.write(tab0, out_ascii, format='fixed_width_two_line')

 
        tab0.write(out_fits,format = 'fits')

    #Plots
    #name = 'Grid_temperature_vs_R23.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)

    fig1, ax1 = plt.subplots()
    ax1.scatter(T_e, R23_composite, marker = '.')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('R_23')
    ax1.set_title('Temperatures_vs_R23')
    ax1.set_xlim(1000,21500)
    pdf_pages.savefig()
     
    fig2, ax2 = plt.subplots()
    ax2.scatter(T_e, O32_composite, marker = '.')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('O_32')
    ax2.set_title('Temperatures_vs_O32')
    ax2.set_xlim(1000,21500)
    pdf_pages.savefig()

    fig3, ax3 = plt.subplots()
    ax3.scatter(R23_composite, com_O_log, marker = '.')
    ax3.set_xlabel('R23')
    ax3.set_ylabel('12+log(O/H) Te')
    ax3.set_title('R23 vs. Composite Metalicity')
    #ax2.set_xlim(1000,21500)
    pdf_pages.savefig()

    fig4, ax4 = plt.subplots()
    ax4.scatter(O32_composite, com_O_log, marker = '.')
    ax4.set_xlabel('O32')
    ax4.set_ylabel('12+log(O/H) Te')
    ax4.set_title('O32 vs. Composite Metalicity')
    #ax2.set_xlim(1000,21500)
    pdf_pages.savefig()
    

    pdf_pages.close()

    

    '''#Histogram
    name = 'Temperature_histogram.pdf'
    pdf_pages = PdfPages(fitspath+name)
    plt.hist(valid_T_e, bins =8)
    plt.xlabel('Temperature (K)')
    #plt.set_ylabel('Spectra')
    plt.title('Preliminary Temperatures')
    pdf_pages.savefig()'''

    #3D plots
    fig_3d= plt.figure(figsize=(10,8))
    ax_3d = plt.axes(projection='3d')
    ax_3d.set_xlabel('R_23')
    ax_3d.set_ylabel('O_32')
    ax_3d.set_zlabel('Temperature')
    ax_3d.set_zlim(4000,26000)
    ax_3d.scatter(R23_composite, O32_composite, T_e, marker='.', linewidths = None)
    plt.show()
    
   

#log(O/H)
#error propagation
'''we have a flux and an sigma flux = flux/signa to noise
propagation distrubution function... monticarlo'''

#3D plots


#Plot O/H values and try the linear plots on the voronoi 14 values 
