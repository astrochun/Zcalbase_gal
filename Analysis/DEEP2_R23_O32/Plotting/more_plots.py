
###PLOTTING RESULTS FROM ANALYSIS
###CALLED EVERYTIME GENERAL FUNCTIONS ARE RUN


###Plotting functions for Zcalbase_gal###


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

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import zoom_and_gauss_general

fitspath_ini='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0502/'
#dataset = 'Double Bin'

lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  

line_type = ['Oxy2', 'Single','Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  


'''
asc_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/Double_Bin_combined_flux_table.tbl'
asc_tab = asc.read(asc_table)

temp_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/Double_Bin_temperatures_metalicity.tbl'
temp_tab = asc.read(temp_table)

verif_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/Double_Bin_verification_master.tbl'
ver_tab = asc.read(verif_table)

R23 = asc_tab['R_23_Average'].data
O32 = asc_tab['O_32_Average'].data
T_e = ver_tab['Temperature'].data
com_O = ver_tab['com_O_log'].data
ID = ver_tab['ID'].data
detect = ver_tab['Detection'].data'''

def ew_plot_R23(fitspath, asc_table, temp_table, verif_table):
    name = fitspath +'equivalent_vs_R23_plots.pdf'
    pdf_pages = PdfPages(name)

    asc_tab = asc.read(asc_table)
    temp_tab = asc.read(temp_table)
    ver_tab = asc.read(verif_table)

    R23 = asc_tab['R_23_Average'].data
    O32 = asc_tab['O_32_Average'].data
    T_e = temp_tab['Temperature'].data
    com_O = temp_tab['com_O_log'].data
    ID = temp_tab['ID'].data
    detect = ver_tab['Detection'].data



    for oo in range(len(lambda0)):
        if line_type[oo] == 'Balmer':
            equ_w = asc_tab['EW_'+str(np.int(lambda0[oo]))+'_abs'].data
            #print equ_w
            #print R23


            fig,ax = plt.subplots()
            ax.scatter(equ_w, R23, marker = '.')
            ax.set_xlabel('Equivalent Width')
            ax.set_ylabel('R23')
            ax.set_title('EW vs. R23   '+str(np.int(lambda0[oo])))
            fig.set_size_inches(8,8)
            fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()



def ew_plot_O32(fitspath, asc_table, temp_table, verif_table):
    name = fitspath +'equivalent_vs_O32_plots.pdf'
    pdf_pages = PdfPages(name)

    
    asc_tab = asc.read(asc_table)
    temp_tab = asc.read(temp_table)
    ver_tab = asc.read(verif_table)

    R23 = asc_tab['R_23_Average'].data
    O32 = asc_tab['O_32_Average'].data
    T_e = temp_tab['Temperature'].data
    com_O = temp_tab['com_O_log'].data
    ID = temp_tab['ID'].data
    detect = ver_tab['Detection'].data

    for oo in range(len(lambda0)):
        if line_type[oo] == 'Balmer':
            equ_w = asc_tab['EW_'+str(np.int(lambda0[oo]))+'_abs'].data
            #print equ_w
            #print R23


            fig,ax = plt.subplots()
            ax.scatter(equ_w, O32, marker= '.')
            ax.set_xlabel('Equivalent Width')
            ax.set_ylabel('O32')
            ax.set_title('EW vs. O32   '+str(np.int(lambda0[oo])))
            fig.set_size_inches(8,8)
            fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()

def R23_vs_O32_color(fitspath, asc_table, temp_table, verif_table):
    name = fitspath +'R23_vs_O32_colormapping.pdf'
    pdf_pages = PdfPages(name)

    
    asc_tab = asc.read(asc_table)
    temp_tab = asc.read(temp_table)
    ver_tab = asc.read(verif_table)

    R23 = asc_tab['R_23_Average'].data
    O32 = asc_tab['O_32_Average'].data
    T_e = temp_tab['Temperature'].data
    com_O = temp_tab['com_O_log'].data
    ID = temp_tab['ID'].data
    detect = ver_tab['Detection'].data

    cm= plt.cm.get_cmap('Blues')
    edge_det = np.where((detect ==1))[0]
    edge_nan = np.where((detect ==0))[0]

    fig1,ax1 = plt.subplots()
    p1= ax1.scatter(R23[edge_det], O32[edge_det], marker= 'o', c=T_e[edge_det], cmap=cm)    #edgecolors = edgecolor
    ax1.scatter(R23[edge_nan], O32[edge_nan], marker= '^')   #, c=T_e, cmap=cm)
    cb = fig1.colorbar(p1)
    cb.set_label('Temperature')
    for aa in range(len(ID)):
        ax1.annotate(ID[aa], (R23[aa], O32[aa]), fontsize = '6')
    ax1.set_xlabel('R23')
    ax1.set_ylabel('O32')
    ax1.set_title('R23 vs. O32 Colormap= Temperature')
    fig1.set_size_inches(8,8)
    fig1.savefig(pdf_pages, format='pdf')

    
    fig2,ax2 = plt.subplots()
    p2 = ax2.scatter(R23[edge_det], O32[edge_det], marker= 'o', c=com_O[edge_det])  #, cmap =cm)#edgecolors = edgecolor
    ax2.scatter(R23[edge_nan], O32[edge_nan], marker= '^')   #, c=com_O, cmap =cm) 
    cb= fig2.colorbar(p2)
    cb.set_label('Metallicity')
    for bb in range(len(ID)):
        ax2.annotate(ID[bb], (R23[bb], O32[bb]), fontsize = '6')
    ax2.set_xlabel('R23')
    ax2.set_ylabel('O32')
    ax2.set_title('R23 vs. O32  Colormap=Metallicity')
    fig2.set_size_inches(8,8)
    fig2.savefig(pdf_pages, format='pdf')





    pdf_pages.close()
    


def hist_for_bin(dataset,asc_table_det3):
    #asc_table = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    asc_tab = asc.read(asc_table_det3)
    name = dataset +'_histograms.pdf'

    pdf_pages = PdfPages(fitspath+name)

    R23 = asc_tab['R23']
    O32 = asc_tab['O32']
    N_bin = asc_tab['N_bin']
    

    for ii in range(np.max(N_bin)+1):
        bin_idx = np.where((N_bin == ii))[0]
        fig, ax = plt.subplots()
        ax.hist(R23[bin_idx])
        ax.set_title('R23 Histogram for Each Bin: Bin'+str(ii))
        pdf_pages.savefig()

        fig.clear()

    pdf_pages.close()
        

    pdf_pages2 = PdfPages(fitspath+'O32'+name)
    for ii in range(np.max(N_bin)+1):
        bin_idx = np.where((N_bin == ii))[0]
        fig, ax = plt.subplots()
        ax.hist(O32[bin_idx])
        ax.set_title('O32 Histogram for Each Bin: Bin'+str(ii))
        pdf_pages2.savefig()

        fig.clear()

    pdf_pages2.close()



def dust_att_plot(combine_flux):
    #name = fitspath+'dust_att_plots.pdf'
    #pdf_pages = PdfPages(name)
    com_asc = asc.read(combine_flux)
    H_gamma_obs = com_asc['Hgamma_Flux_Observed']
    H_beta_obs = com_asc['OIII_4958_Flux_Observed']
    R23 = com_asc['R_23_Average']
    O32 = com_asc['O_32_Average']

    Gambet = H_gamma_obs/H_beta_obs

    fig,ax = plt.subplots()
    ax.scatter(Gambet, O32, marker= '.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('O32')
    ax.set_title('H_gamma/H_beta vs. O32')
    plt.show()
    #fig.set_size_inches(8,8)
    #fig.savefig(pdf_pages, format='pdf')

    fig,ax = plt.subplots()
    ax.scatter(Gambet, R23, marker= '.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('R23')
    ax.set_title('H_gamma/H_beta vs. R23')
    plt.show()
    #fig.set_size_inches(8,8)
    #fig.savefig(pdf_pages, format='pdf')



    #pdf_pages.close()

    
def plotting_individual_for_stacking_image():
    name = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_plots_for_stacking_image.pdf'
    pdf_pages = PdfPages(name)
    
    RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave= header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    scalefactor = 1e-17
    image2d = image2DM/scalefactor
    for ii in range(100,120):
        fig,ax = plt.subplots()
        plt.plot(wave, image2d[ii,:], linewidth = 0.5)
        plt.xlabel('Wavelength')
        plt.ylabel('Intensity Scale factor = 1e-17')
        ax.set_xlim(4100,4500)
                 
        pdf_pages.savefig()
        fig.clear()
    pdf_pages.close()


def plotting_gaussian_curves():
    
    fig, (ax1,ax2,ax3) = plt.subplots(1,3)
    x = np.arange(1,100)
    xbar = 50.0
    s = 15.0
    a = 20.0
    c = 0.0

    singlecurve = np.zeros(len(x))
    singlecurve = zoom_and_gauss_general.gauss(x,xbar,s,a,c)

    #plt.plot(x,singlecurve)
    

    ###Balmer Emission Lines 
    x = np.arange(1,100)
    xbar = 50.0
    s1 = 15.0
    a1 = 20.0
    s2 = 25.0
    a2 = -2.0
    c = 0.0


    doublecurve = np.zeros(len(x))
    doublecurve = zoom_and_gauss_general.double_gauss(x,xbar,s1,a1,c,s2,a2)

    positive = zoom_and_gauss_general.gauss(x,xbar,s1,a1,doublecurve[0])
    negative = zoom_and_gauss_general.gauss(x,xbar,s2,a2,c)
    
    '''plt.plot(x,doublecurve)
    plt.plot(x,positive, linestyle ='--')
    plt.plot(x,negative, linestyle ='--')
    plt.show()'''
    

    ###Oxygen Two Line
    x = np.arange(1,100)
    xbar = 40.0
    s1 = 8.0
    a1 = 20.0
    s2 = 8.0
    a2 = 30.0
    c = 0.0


    oxycurve = np.zeros(len(x))
    oxycurve = oxy2_gauss(x,xbar,s1,a1,c,s2,a2)

    xbar3 = 40.0
    xbar4 = 63.5
    s3 = 8.0
    a3 = 20.0

    s4 = 8.0
    a4 = 30.0
   
    
    positive1 = zoom_and_gauss_general.gauss(x,xbar3,s3,a3,oxycurve[0])
    positive2 = zoom_and_gauss_general.gauss(x,xbar4,s4,a4,oxycurve[0])
    
    '''plt.plot(x,oxycurve)
    plt.plot(x,positive1, 'g', linestyle ='--')
    plt.plot(x,positive2, 'r', linestyle ='--', )
    plt.show()'''

    ax1.plot(x,singlecurve)
    ax2.plot(x,doublecurve)
    ax2.plot(x,positive,'r' ,linestyle ='--')
    ax2.plot(x,negative,'g' ,linestyle ='--')
    ax3.plot(x,oxycurve)
    ax3.plot(x,positive1, 'r', linestyle ='--')
    ax3.plot(x,positive2, 'r', linestyle ='--', )
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax1.set_ylim(-3,25.5)
    ax2.set_ylim(-3,20.5)
    ax3.set_ylim(-3,30.5)
    ax3.set_yticklabels([])
    ax1.set_title('Single Gaussian Curve')
    ax2.set_title('Balmer Fitting with Gaussian Curves')
    ax3.set_title('[OII] Fitting with Gaussian Curves')
    txt1 = '(A)'
    txt2 = '(B)'
    txt3 = '(C)'
    ax1.annotate(txt1, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '10')
    ax2.annotate(txt2, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '10')
    ax3.annotate(txt3, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '10')



    plt.show()

def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    #con1 = 3728.91/3726.16
    con1 = 72.0/45.0
    return a1*np.exp(-(x-xbar)**2/(2*s1**2)) + c + a2*np.exp(-(x-(xbar*con1))**2/(2*s2**2)) 