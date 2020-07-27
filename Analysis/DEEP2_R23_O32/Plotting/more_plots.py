
"""
PLOTTING RESULTS FROM ANALYSIS
Some were written before MSC and are now in MSC

Functions:
-Equivalent Width for R23
-Equivalent Width for O32
-R23 vs O32 Color Map
-Histogram plots for each bin
-Dust attenuation plots (old version before implementing MSC)
-Plotting Individual Spectra over a Wavelength range for Balmer line viewing
-Gaussian Curve Plotting
-Calculating Oxygen II Gaussian
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join
from Zcalbase_gal.Analysis.DEEP2_R23_O32 import zoom_and_gauss_general

from Metallicity_Stack_Commons import lambda0, line_type, line_name

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'


def ew_plot_R23(fitspath, asc_table):
    pdf_pages = PdfPages(join(fitspath, 'equivalent_vs_R23_plots.pdf'))

    asc_tab = asc.read(asc_table)
    R23 = asc_tab['logR23_avg'].data

    for oo in range(len(lambda0)):
        if line_type[oo] == 'Balmer':
            equ_w = asc_tab['EW_'+str(np.int(lambda0[oo]))+'_abs'].data

            fig, ax = plt.subplots()
            ax.scatter(equ_w, R23, marker='.')
            ax.set_xlabel('Equivalent Width')
            ax.set_ylabel('R23')
            ax.set_title('EW vs. R23   '+str(np.int(lambda0[oo])))
            fig.set_size_inches(8,8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()
    pdf_pages.close()


def ew_plot_O32(fitspath, asc_table):
    pdf_pages = PdfPages(join(fitspath,'equivalent_vs_O32_plots.pdf'))

    asc_tab = asc.read(asc_table)
    O32 = asc_tab['logO32_avg'].data

    for oo in range(len(lambda0)):
        if line_type[oo] == 'Balmer':
            equ_w = asc_tab['EW_'+str(np.int(lambda0[oo]))+'_abs'].data

            fig, ax = plt.subplots()
            ax.scatter(equ_w, O32, marker='.')
            ax.set_xlabel('Equivalent Width')
            ax.set_ylabel('O32')
            ax.set_title('EW vs. O32   '+str(np.int(lambda0[oo])))
            fig.set_size_inches(8,8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()
    pdf_pages.close()


def R23_vs_O32_color(fitspath, asc_table, temp_table, verif_table):
    """
    Purpose
    ----------
    Plotting function for R23 and O32 color mapping plots

    Parameters
    ----------
    asc_table   -> combine_flux_ascii
    temp_table  -> derived_properties
    verif_table -> bin_validation_revised
    """

    pdf_pages = PdfPages(join(fitspath, 'R23_vs_O32_colormapping.pdf'))

    asc_tab = asc.read(asc_table)
    temp_tab = asc.read(temp_table)
    ver_tab = asc.read(verif_table)

    R23 = asc_tab['logR23_avg'].data
    O32 = asc_tab['logO32_avg'].data
    T_e = temp_tab['T_e'].data
    com_O = temp_tab['12+log(O/H)'].data
    ID = temp_tab['bin_ID'].data
    detect = ver_tab['Detection'].data

    cm= plt.cm.get_cmap('Blues')
    edge_det = np.where(detect == 1)[0]
    edge_rlimit = np.where(detect == 0.5)[0]

    fig1, ax1 = plt.subplots()
    p1 = ax1.scatter(R23[edge_det], O32[edge_det], marker='o', c=T_e[edge_det], cmap=cm)
    ax1.scatter(R23[edge_rlimit], O32[edge_rlimit], marker='^', c=T_e[edge_rlimit], cmap=cm)
    cb = fig1.colorbar(p1)
    cb.set_label('Temperature')
    for aa in range(len(edge_det)):
        ax1.annotate(ID[edge_det][aa], (R23[edge_det][aa], O32[edge_det][aa]), fontsize='6')
    ax1.set_xlabel('R23')
    ax1.set_ylabel('O32')
    ax1.set_title('R23 vs. O32 Colormap= Temperature')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pdf_pages, format='pdf')
    fig1.clear()

    fig2, ax2 = plt.subplots()
    p2 = ax2.scatter(R23[edge_det], O32[edge_det], marker='o', c=com_O[edge_det])
    ax2.scatter(R23[edge_rlimit], O32[edge_rlimit], marker='^', c=com_O[edge_rlimit])
    cb = fig2.colorbar(p2)
    cb.set_label('Metallicity')
    for bb in range(len(ID)):
        ax2.annotate(ID[bb], (R23[bb], O32[bb]), fontsize='6')
    ax2.set_xlabel('R23')
    ax2.set_ylabel('O32')
    ax2.set_title('R23 vs. O32  Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pdf_pages, format='pdf')
    fig2.clear()

    pdf_pages.close()
    

def hist_for_bin(fitspath, dataset, asc_table_det3):
    asc_tab = asc.read(asc_table_det3)

    pdf_pages = PdfPages(join(fitspath, dataset, '_histograms.pdf')

    R23 = asc_tab['logR23_avg'].data
    O32 = asc_tab['logO32_avg'].data
    N_bin = asc_tab['bin_ID'].data

    number_of_bins = np.max(N_bin) + 1

    for ii in range(number_of_bins):
        bin_idx = np.where(N_bin == ii)[0]
        fig, ax = plt.subplots()
        ax.hist(R23[bin_idx])
        ax.set_title('R23 Histogram for Each Bin: Bin'+str(ii))
        pdf_pages.savefig()

        fig.clear()

    pdf_pages.close()
        

    pdf_pages2 = PdfPages(join(fitspath, 'O32', dataset, '_histograms.pdf')
    for ii in range(number_of_bins):
        bin_idx = np.where(N_bin == ii)[0]
        fig, ax = plt.subplots()
        ax.hist(O32[bin_idx])
        ax.set_title('O32 Histogram for Each Bin: Bin'+str(ii))
        pdf_pages2.savefig()

        fig.clear()

    pdf_pages2.close()


def dust_att_plot(fitspath, combine_flux):
    pdf_pages = PdfPages(join(fitspath, 'dust_attenuation_plots.pdf'))
    com_asc = asc.read(combine_flux)
    H_gamma_obs = com_asc['Hgamma_Flux_Observed']
    H_beta_obs = com_asc['OIII_4958_Flux_Observed']
    R23 = com_asc['R_23_Average']
    O32 = com_asc['O_32_Average']

    Gambet = H_gamma_obs/H_beta_obs

    fig, ax = plt.subplots()
    ax.scatter(Gambet, O32, marker='.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('O32')
    ax.set_title('H_gamma/H_beta vs. O32')
    fig.set_size_inches(8, 8)
    fig.savefig(pdf_pages, format='pdf')

    fig, ax = plt.subplots()
    ax.scatter(Gambet, R23, marker='.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('R23')
    ax.set_title('H_gamma/H_beta vs. R23')
    fig.set_size_inches(8, 8)
    fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()

    
def plotting_individual_for_stacking_image(RestframeMaster, pdf_name, stack_spectra=False):
    if not stack_spectra:
        # name = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_plots_for_stacking_image.pdf'
        # RestframeMaster = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
        spec_range = range(100, 120)
        y_lim = (-2, 2)
        y_text = 1.95
        left = 0.13
    else:
        # RestframeMaster = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/Stacking_Masked_MasterGrid_n_Bins.fits'
        # name = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/composite_plots_data_viz.pdf'
        spec_range = range(27)
        y_lim = (0, 1)
        y_text = 0.95
        left = 0.1

    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    pdf_pages = PdfPages(pdf_name)

    txt0 = r'Intensity ($10^{-17}~{\rm erg}~{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)'

    scalefactor = 1e-17
    image2d = image2DM/scalefactor
    for ii in spec_range:
        fig, ax = plt.subplots()
        plt.plot(wave, image2d[ii, :], linewidth=0.5)
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel(txt0)  # ergs per second per cm^2 per Angstrom
        ax.set_xlim(4250, 4450)
        ax.set_ylim(y_lim)
        ax.axvline(x=4363.21, linewidth=1.0, color='r', linestyle=':')
        ax.axvline(x=4340.46, linewidth=1.0, color='r', linestyle=':')
        ax.text(4363.21, y_text, r'[OIII]$\lambda$4363', va='top', ha='center', rotation=90)
        ax.text(4340.544, y_text, r'H$\gamma$', va='top', ha='center', rotation=90)
        fig.set_size_inches(6, 6)
        plt.subplots_adjust(left=left, right=0.97, bottom=0.1, top=0.97)
        pdf_pages.savefig()
        fig.clear()
    pdf_pages.close()


def plotting_gaussian_curves():
    
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3)
    x = np.arange(1,100)
    xbar = 50.0
    s = 15.0
    a = 20.0
    c = 0.0

    singlecurve = zoom_and_gauss_general.gauss(x, xbar, s, a, c)

    # Balmer Emission Lines
    x = np.arange(1,100)
    xbar = 50.0
    s1 = 15.0
    a1 = 20.0
    s2 = 25.0
    a2 = -2.0
    c = 0.0

    doublecurve = zoom_and_gauss_general.double_gauss(x,xbar,s1,a1,c,s2,a2)

    positive = zoom_and_gauss_general.gauss(x,xbar,s1,a1,doublecurve[0])
    negative = zoom_and_gauss_general.gauss(x,xbar,s2,a2,c)


    # Oxygen Two Line
    x = np.arange(1,100)
    xbar = 40.0
    s1 = 8.0
    a1 = 20.0
    s2 = 8.0
    a2 = 30.0
    c = 0.0

    oxycurve = oxy2_gauss(x,xbar,s1,a1,c,s2,a2)

    xbar3 = 40.0
    xbar4 = 63.5
    s3 = 8.0
    a3 = 20.0

    s4 = 8.0
    a4 = 30.0
    
    positive1 = zoom_and_gauss_general.gauss(x,xbar3,s3,a3,oxycurve[0])
    positive2 = zoom_and_gauss_general.gauss(x,xbar4,s4,a4,oxycurve[0])

    ax1.plot(x, singlecurve)
    ax2.plot(x, doublecurve)
    ax2.plot(x, positive, 'r', linestyle='--')
    ax2.plot(x, negative, 'g', linestyle='--')
    ax3.plot(x, oxycurve)
    ax3.plot(x, positive1, 'r', linestyle='--')
    ax3.plot(x, positive2, 'r', linestyle='--', )
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])
    ax1.set_ylim(-3, 25.5)
    ax2.set_ylim(-3, 20.5)
    ax3.set_ylim(-3, 30.5)
    ax3.set_yticklabels([])
    ax1.set_title('Single Gaussian Curve')
    ax2.set_title('Balmer Fitting with Gaussian Curves')
    ax3.set_title('[OII] Fitting with Gaussian Curves')
    txt1 = '(A)'
    txt2 = '(B)'
    txt3 = '(C)'
    ax1.annotate(txt1, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize='10')
    ax2.annotate(txt2, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize='10')
    ax3.annotate(txt3, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize='10')

    plt.show()


def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    con1 = 72.0/45.0
    return a1*np.exp(-(x-xbar)**2/(2*s1**2)) + c + a2*np.exp(-(x-(xbar*con1))**2/(2*s2**2)) 
