
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from Zcalbase_gal.analysis.deep2_r23_o32 import zoom_and_gauss_general


def r23_vs_o32_color(fitspath, ascii_file, bin_derived_prop_file, verify_file):
    """
    Plotting function for R23 and O32 color mapping plots

    :param fitspath: str. path where files are called from and saved to
    :param ascii_file: str. name of table combine_flux_ascii
    :param bin_derived_prop_file: str. name of table with temperature and
                                  metallicity measurements derived_properties
    :param verify_file: str. name of verification table bin_validation_revised

    PDF File: fitspath + 'R23_vs_O32_colormapping.pdf'
    No Returns
    """

    pdf_file = join(fitspath, 'R23_vs_O32_colormapping.pdf')
    pp = PdfPages(pdf_file)

    ascii_tab = asc.read(ascii_file)
    bin_derived_prop_tab = asc.read(bin_derived_prop_file)
    verify_tab = asc.read(verify_file)

    R23 = ascii_tab['logR23_avg'].data
    O32 = ascii_tab['logO32_avg'].data
    T_e = bin_derived_prop_tab['T_e'].data
    com_O = bin_derived_prop_tab['12+log(O/H)'].data
    ID = bin_derived_prop_tab['bin_ID'].data
    detect = verify_tab['Detection'].data

    cm = plt.cm.get_cmap('Blues')
    edge_det = np.where(detect == 1)[0]
    edge_rlimit = np.where(detect == 0.5)[0]

    fig1, ax1 = plt.subplots()
    p1 = ax1.scatter(R23[edge_det], O32[edge_det], marker='o',
                     c=T_e[edge_det], cmap=cm)
    ax1.scatter(R23[edge_rlimit], O32[edge_rlimit], marker='^',
                c=T_e[edge_rlimit], cmap=cm)
    cb = fig1.colorbar(p1)
    cb.set_label('Temperature')
    for aa in range(len(edge_det)):
        ax1.annotate(ID[edge_det][aa], (R23[edge_det][aa], O32[edge_det][aa]),
                     fontsize='6')
    ax1.set_xlabel('R23')
    ax1.set_ylabel('O32')
    ax1.set_title('R23 vs. O32 Colormap= Temperature')
    fig1.set_size_inches(8, 8)
    fig1.savefig(pp, format='pdf')
    fig1.clear()

    fig2, ax2 = plt.subplots()
    p2 = ax2.scatter(R23[edge_det], O32[edge_det], marker='o',
                     c=com_O[edge_det])
    ax2.scatter(R23[edge_rlimit], O32[edge_rlimit], marker='^',
                c=com_O[edge_rlimit])
    cb = fig2.colorbar(p2)
    cb.set_label('Metallicity')
    for bb in range(len(ID)):
        ax2.annotate(ID[bb], (R23[bb], O32[bb]), fontsize='6')
    ax2.set_xlabel('R23')
    ax2.set_ylabel('O32')
    ax2.set_title('R23 vs. O32  Colormap=Metallicity')
    fig2.set_size_inches(8, 8)
    fig2.savefig(pp, format='pdf')
    fig2.clear()
    pp.close()


def hist_for_bin(fitspath, dataset, asc_table_det3):
    """
    Produces a pdf file with plots of histograms to check the distribution
    of individual galaxies in bins based on R23 and O32

    :param fitspath: str. path where files are called from and saved to
    :param dataset: str. keyword that specifies which binning method is used
    :param asc_table_det3: str. ascii table created by binning code
                            (ie. bin_info.tbl)

    PDF File: fitspath + 'bin_histograms.pdf'
    No Returns
    """
    asc_tab = asc.read(asc_table_det3)

    pdf_file1 = join(fitspath, 'bin_histograms.pdf')
    pp = PdfPages(pdf_file1)

    R23 = asc_tab['logR23_avg'].data
    O32 = asc_tab['logO32_avg'].data
    N_bin = asc_tab['bin_ID'].data

    number_of_bins = np.max(N_bin) + 1
    for ii in range(number_of_bins):
        bin_idx = np.where(N_bin == ii)[0]
        fig, ax = plt.subplots()
        ax.hist(R23[bin_idx])
        ax.set_title('R23 Histogram for Each Bin: Bin' + str(ii))
        pp.savefig()

        fig.clear()

    pp.close()

    pdf_file2 = join(fitspath, 'O32_histogram.pdf')
    pp2 = PdfPages(pdf_file2)
    for ii in range(number_of_bins):
        bin_idx = np.where(N_bin == ii)[0]
        fig, ax = plt.subplots()
        ax.hist(O32[bin_idx])
        ax.set_title('O32 Histogram for Each Bin: Bin' + str(ii))
        pp2.savefig()

        fig.clear()

    pp2.close()

    
def plotting_individual_for_stacking_image(RestframeMaster, pdf_file,
                                           stack_spectra=False):
    """
    Produces pdf file of spectra plots for either the binned data or
    individual data.
    Names commented out to refer to in future for what files are required.
    Responsible for the plots for Stacking Spectra Figure (from thesis)

    :param RestframeMaster: str. master grid file
    :param pdf_file: str. name of outputted file
    :param stack_spectra: bool. Default = False determines if plotting
                          individual spectra (False) or stacked spectra (True)

    PDF File: pdf_file
    No returns
    """
    if not stack_spectra:
        # name = fitspath_ini + 'individual_plots_for_stacking_image.pdf'
        # RestframeMaster = fitspath_ini + '/Master_Grid.fits'
        spec_range = range(100, 120)
        y_lim = (-2, 2)
        y_text = 1.95
        left = 0.13
    else:
        # RestframeMaster = fitspath + 'Stacking_Masked_MasterGrid_n_Bins.fits'
        # name = fitspath + 'composite_plots_data_viz.pdf'
        spec_range = range(27)
        y_lim = (0, 1)
        y_text = 0.95
        left = 0.1

    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    pp = PdfPages(pdf_file)

    txt0 = r'Intensity ($10^{-17}~{\rm erg}~' \
           r'{\rm s}^{-1}~{\rm cm}^{-2}~\AA^{-1}$)'
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
        ax.axvline(x=4340.544, linewidth=1.0, color='r', linestyle=':')
        ax.text(4363.21, y_text, r'[OIII]$\lambda$4363', va='top', ha='center',
                rotation=90)
        ax.text(4340.46, y_text, r'H$\gamma$', va='top', ha='center',
                rotation=90)
        fig.set_size_inches(6, 6)
        plt.subplots_adjust(left=left, right=0.97, bottom=0.1, top=0.97)
        pp.savefig()
        fig.clear()
    pp.close()


def plotting_gaussian_curves():
    """
    Plots the single, double, and oxygen gaussian curves used
    to fit the binned spectra.

    No returns
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    x = np.arange(1, 100)
    xbar = 50.0
    s = 15.0
    a = 20.0
    c = 0.0
    singlecurve = zoom_and_gauss_general.gauss(x, xbar, s, a, c)

    # Balmer Emission Lines
    x = np.arange(1, 100)
    xbar = 50.0
    s1 = 15.0
    a1 = 20.0
    s2 = 25.0
    a2 = -2.0
    c = 0.0

    doublecurve = zoom_and_gauss_general.double_gauss(x, xbar,
                                                      s1, a1, c, s2, a2)

    positive = zoom_and_gauss_general.gauss(x, xbar, s1, a1, doublecurve[0])
    negative = zoom_and_gauss_general.gauss(x, xbar, s2, a2, c)

    # Oxygen Two Line
    x = np.arange(1, 100)
    xbar = 40.0
    s1 = 8.0
    a1 = 20.0
    s2 = 8.0
    a2 = 30.0
    c = 0.0

    oxycurve = oxy2_gauss(x, xbar, s1, a1, c, s2, a2)

    xbar3 = 40.0
    xbar4 = 63.5
    s3 = 8.0
    a3 = 20.0

    s4 = 8.0
    a4 = 30.0

    positive1 = zoom_and_gauss_general.gauss(x, xbar3, s3, a3, oxycurve[0])
    positive2 = zoom_and_gauss_general.gauss(x, xbar4, s4, a4, oxycurve[0])

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
    ax1.annotate(txt1, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')
    ax2.annotate(txt2, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')
    ax3.annotate(txt3, [0.95, 0.95], xycoords='axes fraction', va='top',
                 ha='right', fontsize='10')

    plt.show()


def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    """
    Calculates the gaussian curve used to fit the OII emission.

    :param x: array. the emission line that is being fitted
    :param xbar, s1, a1, c, s2, a2 : int. guess parameters for the fit

    Returns Oxygen II gaussian fit
    """

    con1 = 72.0/45.0
    return a1 * np.exp(-(x - xbar)**2/(2 * s1**2)) + c + a2 * \
           np.exp(-(x - (xbar * con1))**2/(2 * s2**2))
