# This code generates a pdf file with all the balmer plots for a given bin

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from chun_codes.cardelli import *


fitspath_ini='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
# fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/'
# dataset = 'Double Bin'

lambda0 = [3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Oxy2', 'Single','Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727', 'NeIII', 'HeI', '3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import zoom_and_gauss_general as zm
# scalefact= zm.zoom_gauss_plot.scalefact
scalefact = 1e-17


def balmer_graphs(fitspath, nrow, ncol, Stack_name, combine_flux_tab, out_pdf):

    # Stack_name= '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Stacking_Masked_MasterGrid_singleDouble_Bin.fits'
 
    stack2D, header = fits.getdata(Stack_name, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    dispersion = header['CDELT1']

    combine_asc = asc.read(combine_flux_tab)

    ID = combine_asc['ID']

    B_xbar = combine_asc['HBETA_X_bar']
    B_sp = combine_asc['HBETA_Pos_Sig']
    B_ap = combine_asc['HBETA_Pos_Amp']
    B_con = combine_asc['HBETA_Const']
    B_sn = combine_asc['HBETA_Neg_Sig']
    B_an = combine_asc['HBETA_Neg_Amp']

    G_xbar = combine_asc['Hgamma_X_bar']
    G_sp = combine_asc['Hgamma_Pos_Sig']
    G_ap = combine_asc['Hgamma_Pos_Amp']
    G_con = combine_asc['Hgamma_Const']
    G_sn = combine_asc['Hgamma_Neg_Sig']
    G_an = combine_asc['Hgamma_Neg_Amp']

    D_xbar = combine_asc['HDELTA_X_bar']
    D_sp = combine_asc['HDELTA_Pos_Sig']
    D_ap = combine_asc['HDELTA_Pos_Amp']
    D_con = combine_asc['HDELTA_Const']
    D_sn = combine_asc['HDELTA_Neg_Sig']
    D_an = combine_asc['HDELTA_Neg_Amp']

    # out_pdf = fitspath + '/AllBalmer_line_fitting.pdf'
    # Still need to incoorporate various working waves and importing xo, ynorm
    pdfpages = PdfPages(out_pdf)
    nrows = 3
    ncols = 3

    for ii in range(len(ID)):
   
        if ii % nrows == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact
        dx = wave[2] - wave[1]

        # Beta
        working_wave_beta = 4861.32
        Bx_sigsnip = np.where(np.abs((wave - working_wave_beta))/B_sp[ii] <= 2.5)[0]
        Bgauss0 = zm.double_gauss(wave, B_xbar[ii], B_sp[ii], B_ap[ii], B_con[ii], B_sn[ii], B_an[ii])
        Bneg0 = zm.gauss(wave, B_xbar[ii], B_sn[ii], B_an[ii], B_con[ii])
        Bgauss0_diff = Bgauss0 - Bneg0
        By_norm_diff = y_norm[Bx_sigsnip] - Bneg0[Bx_sigsnip]

        # Residuals
        Bx_sigsnip_2 = np.where(np.abs((wave - working_wave_beta))/B_sp[ii] <= 3.0)[0]
        Bresid = y_norm[Bx_sigsnip_2] - Bgauss0[Bx_sigsnip_2] + B_con[ii]

        # Fluxes
        Bflux_g = np.sum(Bgauss0_diff*dx)
        Bflux_s = np.sum(By_norm_diff*dx)
        
        # Gamma
        working_wave_gamma = 4340.46
        Gx_sigsnip = np.where(np.abs((wave - working_wave_gamma))/G_sp[ii] <= 2.5)[0]
        Ggauss0 = zm.double_gauss(wave, G_xbar[ii], G_sp[ii], G_ap[ii], G_con[ii], G_sn[ii], G_an[ii])
        Gneg0 = zm.gauss(wave, G_xbar[ii], G_sn[ii], G_an[ii], G_con[ii])
        Ggauss0_diff = Ggauss0 - Gneg0
        Gy_norm_diff = y_norm[Gx_sigsnip] - Gneg0[Gx_sigsnip]

        # Residuals
        Gx_sigsnip_2 = np.where(np.abs((wave - working_wave_gamma))/G_sp[ii] <= 3.0)[0]
        Gresid = y_norm[Gx_sigsnip_2] - Ggauss0[Gx_sigsnip_2] + G_con[ii]

        # Fluxes
        Gflux_g = np.sum(Ggauss0_diff*dx)
        Gflux_s = np.sum(Gy_norm_diff*dx)

        # Delta
        working_wave_delta = 4101.73
        Dx_sigsnip = np.where(np.abs((wave - working_wave_delta))/D_sp[ii] <= 2.5 )[0]
        Dgauss0 = zm.double_gauss(wave, D_xbar[ii], D_sp[ii], D_ap[ii], D_con[ii], D_sn[ii], D_an[ii])
        Dneg0 = zm.gauss(wave, D_xbar[ii], D_sn[ii], D_an[ii], D_con[ii])
        Dgauss0_diff = Dgauss0 - Dneg0
        Dy_norm_diff = y_norm[Dx_sigsnip] - Dneg0[Dx_sigsnip]

        # Residuals
        Dx_sigsnip_2 = np.where(np.abs((wave - working_wave_delta))/D_sp[ii] <= 3.0)[0]
        Dresid = y_norm[Dx_sigsnip_2] - Dgauss0[Dx_sigsnip_2] + D_con[ii]

        # Fluxes
        Dflux_g = np.sum(Dgauss0_diff*dx)
        Dflux_s = np.sum(Dy_norm_diff*dx)

        row = ii % nrows

        txt0 = r'ID: %i' % (ID[ii]) + '\n'
        txt0 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (B_sp[ii], B_sn[ii]) + '\n'
        txt0 += 'F_G: %.3f F_S: %.3f' % (Bflux_g, Bflux_s)
        # txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) +
    
        ax_arr[row][2].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][2].plot(wave, Bgauss0, 'm', linewidth=0.25, label='Beta Fit')
        ax_arr[row][2].set_xlim(working_wave_beta-50, working_wave_beta+50)
        # ax_arr[row][2].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2,fontsize = 3)
        ax_arr[row][2].annotate(txt0, [0.95, 0.95], xycoords='axes fraction', va='top', ha='right', fontsize='5')
        ax_arr[row][2].plot(wave[Bx_sigsnip_2], Bresid, 'r', linestyle='dashed', linewidth=0.2, label='Residuals')

        txt1 = r'ID: %i' % (ID[ii]) + '\n'
        txt1 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (G_sp[ii], G_sn[ii]) + '\n'
        txt1 += 'F_G: %.3f F_S: %.3f' % (Gflux_g, Gflux_s)
    
        ax_arr[row][1].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][1].plot(wave, Ggauss0, 'm', linewidth=0.25, label='Gamma Fit')
        ax_arr[row][1].set_xlim(working_wave_gamma - 50, working_wave_gamma + 50)
        # ax_arr[row][1].legend(bbox_to_anchor=(0.25, 0.1), borderaxespad=0, ncol=2, fontsize=3)
        ax_arr[row][1].annotate(txt1, [0.95, 0.95], xycoords='axes fraction', va='top', ha='right', fontsize='5')
        ax_arr[row][1].plot(wave[Gx_sigsnip_2], Gresid, 'r', linestyle='dashed', linewidth=0.2, label='Residuals')

        txt2 = r'ID: %i' % (ID[ii]) + '\n'
        txt2 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (D_sp[ii], D_sn[ii]) + '\n'
        txt2 += 'F_G: %.3f F_S: %.3f' % (Dflux_g, Dflux_s)

        ax_arr[row][0].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][0].plot(wave, Dgauss0, 'm', linewidth=0.25, label='Detla Fit')
        ax_arr[row][0].set_xlim(working_wave_delta - 50, working_wave_delta + 50)

        ax_arr[row][0].set_ylim(0, 1.5)
        
        # ax_arr[row][0].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize=3)
        ax_arr[row][0].annotate(txt0, [0.95, 0.95], xycoords='axes fraction', va='top', ha='right', fontsize='5')
        ax_arr[row][0].plot(wave[Dx_sigsnip_2], Dresid, 'r', linestyle='dashed', linewidth=0.2, label='Residuals')

        ax_arr[row][0].set_yticklabels([0, 0.5, 1, 1.5])
        ax_arr[row][1].set_yticklabels([])  # sets y-tick labels
        ax_arr[row][2].set_yticklabels([])

        if row != nrows - 1 and ii != stack2D.shape[0] - 1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])

        if ii % nrows == nrows-1:
            fig.savefig(pdfpages, format='pdf')

    pdfpages.close()


'''
    # H_Beta_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/HBETA_Balmer_fitting_values.tbl'
    H_Beta = asc.read(HBETA_table)
    # H_Gamma_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Hgamma_Balmer_fitting_values.tbl'
    H_Gamma = asc.read(HGAMMA_table)
    # H_Delta_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/HDELTA_Balmer_fitting_values.tbl'
    H_Delta = asc.read(HDELTA_table)

    ID = H_Beta['ID']

    B_xbar = H_Beta['X_bar']
    B_sp = H_Beta['Pos_Sig']
    B_ap = H_Beta['Pos_Amp']
    B_con = H_Beta['Const']
    B_sn = H_Beta['Neg_Sig']
    B_an = H_Beta['Neg_Amp']

    G_xbar = H_Gamma['X_bar']
    G_sp = H_Gamma['Pos_Sig']
    G_ap = H_Gamma['Pos_Amp']
    G_con = H_Gamma['Const']
    G_sn = H_Gamma['Neg_Sig']
    G_an = H_Gamma['Neg_Amp']

    D_xbar = H_Delta['X_bar']
    D_sp = H_Delta['Pos_Sig']
    D_ap = H_Delta['Pos_Amp']
    D_con = H_Delta['Const']
    D_sn = H_Delta['Neg_Sig']
    D_an = H_Delta['Neg_Amp']'''
