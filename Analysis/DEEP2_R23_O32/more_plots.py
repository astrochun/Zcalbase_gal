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


fitspath_ini='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/'
dataset = 'Double Bin'

lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  

line_type = ['Oxy2', 'Single','Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  



asc_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/Double_Bin_combined_flux_table.tbl'
asc_tab = asc.read(asc_table)

temp_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_0206/Double_Bin_temperatures_metalicity.tbl'
temp_tab = asc.read(temp_table)
R23 = asc_tab['R_23_Average'].data
O32 = asc_tab['O_32_Average'].data
T_e = temp_tab['Temperature'].data
com_O = temp_tab['com_O_log'].data
ID = temp_tab['ID'].data

def ew_plot_R23():
    name = fitspath +'equivalent_vs_R23_plots.pdf'
    pdf_pages = PdfPages(name)
    



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



def ew_plot_O32():
    name = fitspath +'equivalent_vs_O32_plots.pdf'
    pdf_pages = PdfPages(name)

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

def R23_vs_O32_color():
    name = fitspath +'R23_vs_O32_colormapping.pdf'
    pdf_pages = PdfPages(name)

    cm= plt.cm.get_cmap('Blues')
    
    fig1,ax1 = plt.subplots()
    p1= ax1.scatter(R23, O32, marker= 'o', c=T_e, cmap=cm)
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
    p2 = ax2.scatter(R23, O32, marker= 'o', c=com_O, cmap =cm)
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
    

