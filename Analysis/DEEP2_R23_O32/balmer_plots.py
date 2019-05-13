###This code generates a pdf file with all the balmer plots for a given bin###


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
fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/'
#dataset = 'Double Bin'

lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  

line_type = ['Oxy2', 'Single','Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import zoom_and_gauss_general as zm
#scalefact= zm.zoom_gauss_plot.scalefact
scalefact = 1e-17


Stack_name= '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Stacking_Masked_MasterGrid_singleDouble_Bin.fits'
 
stack2D, header = fits.getdata(Stack_name,header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
dispersion = header['CDELT1']




H_Beta_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/HBETA_Balmer_fitting_values.tbl'
H_Beta = asc.read(H_Beta_tab)
H_Gamma_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Hgamma_Balmer_fitting_values.tbl'
H_Gamma = asc.read(H_Gamma_tab)
H_Delta_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/HDELTA_Balmer_fitting_values.tbl'
H_Delta = asc.read(H_Delta_tab)

ID = H_Beta['ID']

B_xbar = H_Beta['X_bar']
B_sp= H_Beta['Pos_Sig']
B_ap= H_Beta['Pos_Amp']
B_con= H_Beta['Const']
B_sn= H_Beta['Neg_Sig']
B_an= H_Beta['Neg_Amp']

G_xbar = H_Gamma['X_bar']
G_sp= H_Gamma['Pos_Sig']
G_ap= H_Gamma['Pos_Amp']
G_con= H_Gamma['Const']
G_sn= H_Gamma['Neg_Sig']
G_an= H_Gamma['Neg_Amp']

D_xbar = H_Delta['X_bar']
D_sp= H_Delta['Pos_Sig']
D_ap= H_Delta['Pos_Amp']
D_con= H_Delta['Const']
D_sn= H_Delta['Neg_Sig']
D_an= H_Delta['Neg_Amp']





###Still need to incoorporate various working waves and importing xo, ynorm
pdfpages = PdfPages(fitspath + '/AllBalmer_line_fitting.pdf')
nrows = 3
ncols = 3

for ii in range(len(ID)):
   
    if ii % nrows ==0: fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)

    y0 = stack2D[ii]
    y_norm = y0/scalefact


    
    ##Beta
    working_wave_beta = 4861.32
    Bx_sigsnip = np.where(np.abs((wave-working_wave_beta))/B_sp[ii]<=2.5 )[0]
    Bgauss0 = zm.double_gauss(wave, B_xbar[ii], B_sp[ii], B_ap[ii], B_con[ii], B_sn[ii], B_an[ii])
    Bneg0   = zm.gauss(wave, B_xbar[ii], B_sn[ii], B_an[ii], B_con[ii])
    Bgauss0_diff = Bgauss0 - Bneg0
    By_norm_diff = y_norm[Bx_sigsnip]-Bneg0[Bx_sigsnip]

    ##Gamma
    working_wave_gamma = 4340.46
    Gx_sigsnip = np.where(np.abs((wave-working_wave_gamma))/G_sp[ii]<=2.5 )[0]
    Ggauss0 = zm.double_gauss(wave, G_xbar[ii], G_sp[ii], G_ap[ii], G_con[ii], G_sn[ii], G_an[ii])
    Gneg0   = zm.gauss(wave, G_xbar[ii], G_sn[ii], G_an[ii], G_con[ii])
    Ggauss0_diff = Ggauss0 - Gneg0
    Gy_norm_diff = y_norm[Gx_sigsnip]-Gneg0[Gx_sigsnip]


    ##Delta
    working_wave_delta =  4101.73
    Dx_sigsnip = np.where(np.abs((wave-working_wave_delta))/D_sp[ii]<=2.5 )[0]
    Dgauss0 = zm.double_gauss(wave, D_xbar[ii], D_sp[ii], D_ap[ii], D_con[ii], D_sn[ii], D_an[ii])
    Dneg0   = zm.gauss(wave, D_xbar[ii], D_sn[ii], D_an[ii], D_con[ii])
    Dgauss0_diff = Dgauss0 - Dneg0
    Dy_norm_diff = y_norm[Dx_sigsnip]-Dneg0[Dx_sigsnip]

    row = ii % nrows

    
    ax_arr[row][0].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
    ax_arr[row][0].plot(wave,Bgauss0, 'm', linewidth= 0.25, label= 'Beta Fit')
    ax_arr[row][0].set_xlim(working_wave_beta-50, working_wave_beta+50)
    #ax_arr[row][0].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)
    txt0 = r'ID: %i, Beta Emission' % (ID[ii])
    ax_arr[row][0].annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
    
    ax_arr[row][1].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
    ax_arr[row][1].plot(wave,Ggauss0, 'm', linewidth= 0.25, label= 'Gamma Fit')
    ax_arr[row][1].set_xlim(working_wave_gamma-50, working_wave_gamma+50)
    #ax_arr[row][1].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)
    txt0 = r'ID: %i, Gamma Emission' % (ID[ii])
    ax_arr[row][1].annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')

    ax_arr[row][2].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
    ax_arr[row][2].plot(wave,Dgauss0, 'm', linewidth= 0.25, label= 'Detla Fit')
    ax_arr[row][2].set_xlim(working_wave_delta-50, working_wave_delta+50)
    #ax_arr[row][2].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)
    txt0 = r'ID: %i, Delta Emission' % (ID[ii])
    ax_arr[row][2].annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')

    if ii% nrows == nrows-1: fig.savefig(pdfpages, format='pdf')

pdfpages.close()




    ###Will Need to go back and add residual calculation if wanted in plots t_ax.plot(x0[x_sigsnip_2],resid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')
    #t_ax.plot(wave,lineflag,'g', linestyle='dashed', linewidth = 0.1, label = 'Lineflag')
    #t_ax.legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)
    #t_ax.set_xlim(x1+50,x2-50)
