import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from ..log_commons import log_stdout


def color_for_bin(fitspath, bin_info_file, pdf_file):
    """
    Plots for R23 and O32 for Voronoi analysis for each bin
    NOT CALLED IN GENERAL FUNCTIONS

    :param fitspath: str. path where files are called from and saved to
    :param bin_info_file: str. table created by binning code
    :param pdf_file: str. name of pdf file produced

    PDF File: fitspath + pdf_file
    No returns
    """

    pp = PdfPages(join(fitspath, pdf_file))

    bin_info_tab = asc.read(join(fitspath, bin_info_file))
    targetSN = 14
    xBar = bin_info_tab['logR23_avg']
    yBar = bin_info_tab['logO32_avg']
    xnode = bin_info_tab['logR23_min']
    ynode = bin_info_tab['logO32_min']
    area = bin_info_tab['N_stack']

    w = area == 1
    plt.clf()
    plt.subplot(211)
    rnd = np.argsort(np.random.random(xnode.size))  # Randomize bin colors
    # do not rescale after imshow()
    plt.plot(xnode, ynode, '+w', scalex=False, scaley=False)
    plt.xlabel('R (arcsec)')
    plt.ylabel('R (arcsec)')
    plt.title('Map of Voronoi bins')

    plt.subplot(212)
    rad = np.sqrt(xBar**2 + yBar**2)  # Use centroids, NOT generators
    plt.plot(rad[~w], sn[~w], 'or')
    plt.xlabel('R (arcsec)')
    plt.ylabel('Bin S/N')
    plt.axis([np.min(rad), np.max(rad), 0, np.max(sn)])  # x0, x1, y0, y1
    if np.sum(w) > 0:
        plt.plot(rad[w], sn[w], 'xb', label='single spaxels')
    plt.axhline(targetSN)
    plt.legend()
    plt.pause(1)
    pp.savefig()
    pp.close()


def r23vso32_plot(fitspath, bin_info_file, bin_derived_prop_file, pdf_file,
                  log=None):
    """
    Plotting R23 vs O32 with a color map for
    metallicity and then for temperature

    :param fitspath: str. path where files are called from and saved to
    :param bin_info_file: str. table created by binning code
    :param bin_derived_prop_file: str. table holding metallicity and
                                  temperature measurements
    :param pdf_file: str. name of pdf file produced
    :param log: LogClass or logging object

    PDF File: fitspath + pdf_file
    No returns
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    # pp = 'R23vsO32_color_comandavg.pdf'
    pp = PdfPages(join(fitspath, pdf_file))

    bin_info_tab = asc.read(join(fitspath, bin_info_file))

    bin_derived_prop_tab = asc.read(bin_derived_prop_file)

    # logR23 = asc_table['log(R23)']
    # logO32 = asc_table['log(O32)']
    com_O_log = bin_derived_prop_tab['12+log(O/H)']
    R23_com = bin_derived_prop_tab['R23_Composite']
    O32_com = bin_derived_prop_tab['O32_Composite']
    n_gal = bin_derived_prop_tab['N_Stack']
    
    xBar = bin_info_tab['logR23_avg']
    yBar = bin_info_tab['logO32_avg']

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    vmin = np.nanmin(com_O_log)  # temp)
    vmax = np.nanmax(com_O_log[com_O_log != np.inf])   # temp
    c = com_O_log  # temp

    log.info(f"min, max, c: {vmin}, {vmax}, {c}")
    scat = plt.scatter(xBar, yBar, c=c, vmin=vmin, vmax=vmax,
                       marker='.', cmap=cm)
    plt.scatter(R23_com, O32_com, c=c, vmin=vmin, vmax=vmax,
                marker='*', cmap=cm)
    cax = plt.axes([0.9, 0.15, 0.8, 0.75])
    plt.colorbar(scat, cax=cax, orientation='vertical',
                 label='Oxygen Abundance')
    for rr in range(len(n_gal)):
        ax1.annotate(str(n_gal[rr]), (xBar[rr], yBar[rr]), xycoords='data',
                     fontsize=5)
    for oo in range(len(n_gal)):
        ax1.annotate(str(n_gal[oo]), (R23_com[oo], O32_com[oo]),
                     xycoords='data', fontsize=5, fontweight='bold')
    ax1.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax1.set_xlabel(r'log($R_{23}$)')
    ax1.set_ylabel(r'log($O_{32}$)')

    pp.savefig()
    pp.close()

    log.info("finished.")
