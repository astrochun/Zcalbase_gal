import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
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

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

#Voronoi10
#Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi10/Voronoi_combined_flux_table.tbl'
data1 = asc.read(spectral)
'''
#Grid
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_combined_flux_table.tbl'
data1 = asc.read(spectral)'''

#Voronoi14
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_combined_flux_table.tbl'
data1 = asc.read(spectral)


OIII5007 = data1['OIII_5007_Flux_Observed'].data
OIII4959 = data1['OIII_4958_Flux_Observed'].data
OIII4363 = data1['OIII_4363_Flux_Observed'].data
HBETA    = data1['HBETA_Flux_Observed'].data
OII3727  = data1['OII_3727_Flux_Observed'].data
R23      = data1['R_23_Average'].data
O32      = data1['O_32_Average'].data
N_Galaxy = data1['N_Galaxies'].data

SN_5007       = data1['OIII_5007_S/N'].data
SN_4959       = data1['OIII_4958_S/N'].data
SN_4363       = data1['OIII_4363_S/N'].data
SN_HBETA      = data1['HBETA_S/N'].data
SN_3727       = data1['OII_3727_S/N'].data

#Constants

a = 13205
b = 0.92506
c = 0.98062
#x = S/N([O iii]4363)

#Calculations on page 37 of paper

def R_calculation():
  
    R_value = OIII4363/(OIII5007+OIII4959)
    return R_value

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)      #np.zeros(len(OIII5007))
  
    return T_e


def metalicity_calculation(T_e):
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
    '''t = 1e-4 * Te
    x = 1e-4 * 1e3 * t**(-0.5)

    dist_t = 1e-4 * dist_Te
    dist_x = (1e-4 * 1e3 * dist_t**(-0.5))

    lOIIIHb      = np.log10((OIII4959 + OIII5007)/HB)
    dist_lOIIIHb = np.log10((OIII4959_dist + OIII5007_dist)/HB_dist)

    # T2-T3 RELATION
    if andrews == False:
        print '### Using Izotov et al. (2006) T2-T3 relation'
        t2   = -0.577 + t*(2.065-0.498*t)
        mark = np.where(t >= 2.075)[0]
        if len(mark) > 0:
            t2[mark] = 1.57067
    else:
        print '### Using Andrews & Martini (2013) T2-T3 relation'     
        t2   = 0.7*t + (3000.-1300.)/1e4
        mark = np.where(t >= 2.0)
        if len(mark) > 0: 
            t2[mark] = 1.57
    #endelse

    x2 = 1e-4 * 1e3 * t2**(-0.5)
    TOII = t2 * 1e4'''

    return O_s_ion , O_d_ion, com_O_log
def run_function():
    R_value = R_calculation()
    T_e = temp_calculation(R_value)
    O_s_ion, O_d_ion, com_O_log = metalicity_calculation(T_e)
    
    print R_value
    print T_e
    print min(T_e), max(T_e)
    print O_s_ion, O_d_ion

   
    #valid_T_e = [7284.30776602, 7718.16922375, 8605.17485391, 8609.03697828, 10320.26788653, 10332.93036254, 10773.16786121, 11076.54040316, 11579.32008455, 11774.61655556, 11956.34906994,  12273.56242391, 12482.29559535, 15202.11770306, 16070.1413722, 20520.22534753]

    #Ascii Table
    out_ascii = fitspath+ '/Voronoi14_temperatures_metalicity_asc_table.tbl'
    if not exists(out_ascii):
        n=  ('R_23_Average', 'O_32_Average', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959', 'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727', 'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log')
        
        tab0 = Table([ R23, O32, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, T_e, O_s_ion, O_d_ion, com_O_log], names=n)
        asc.write(tab0, out_ascii, format='fixed_width_two_line')
    
    #Plots
    name = 'Voronoi14_temperature_vs_R23.pdf'
    pdf_pages = PdfPages(fitspath+name)

    fig1, ax1 = plt.subplots()
    ax1.scatter(T_e, R23, marker = '.')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('R_23')
    ax1.set_title('Temperatures_vs_R23')
    ax1.set_xlim(1000,21500)
    pdf_pages.savefig()
     
    fig2, ax2 = plt.subplots()
    ax2.scatter(T_e, O32, marker = '.')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('O_32')
    ax2.set_title('Temperatures_vs_O32')
    ax2.set_xlim(1000,21500)
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
    ax_3d.scatter(R23, O32, T_e, marker='.', linewidths = None)
    plt.show()
    
   

#log(O/H)
#error propagation
'''we have a flux and an sigma flux = flux/signa to noise
propagation distrubution function... monticarlo'''

#3D plots
