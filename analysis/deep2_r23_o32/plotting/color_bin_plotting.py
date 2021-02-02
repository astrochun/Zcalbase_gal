import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from ..log_commons import log_stdout


def color_for_bin(fitspath, bin_info, pdf_file):
    """
    Plots for R23 and O32 for Voronoi analysis for each bin
    NOT CALLED IN GENERAL FUNCTIONS

    param: fitspath: str. path where files are called from and saved to
    param: bin_info: str. table created by binning code
    param: pdf_file: str. name of pdf file produced

    PDF File: fitspath + pdf_file
    No returns
    """

    pdf_pages = PdfPages(join(fitspath, pdf_file))
    asc_table = asc.read(join(fitspath, bin_info))
    targetSN = 14
    xBar = asc_table['logR23_avg']
    yBar = asc_table['logO32_avg']
    xnode = asc_table['logR23_min']
    ynode = asc_table['logO32_min']
    area = asc_table['N_stack']

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
    pdf_pages.savefig()
    pdf_pages.close()


def r23vso32_plot(fitspath, bin_info, temp_tab, pdf_name, log=None):
    """
    Plotting R23 vs O32 with a color map for
    metallicity and then for temperature

    :param fitspath: str. path where files are called from and saved to
    :param bin_info: str. table created by binning code
    :param temp_tab: str. table holding metallicity and
                     temperature measurements
    :param pdf_name: str. name of pdf file produced
    :param log: LogClass or logging object

    PDF File: fitspath + pdf_name
    No returns
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    # pdf_name = 'R23vsO32_color_comandavg.pdf'
    pdf_pages = PdfPages(join(fitspath, pdf_name))

    asc_table1 = asc.read(bin_info)
    temp_table = asc.read(temp_tab)

    # logR23 = asc_table['log(R23)']
    # logO32 = asc_table['log(O32)']
    com_O_log = temp_table['12+log(O/H)']
    R23_com = temp_table['R23_Composite']
    O32_com = temp_table['O32_Composite']
    n_gal = temp_table['N_Stack']
    
    xBar = asc_table1['logR23_avg']
    yBar = asc_table1['logO32_avg']

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


    pdf_pages.savefig()
    pdf_pages.close()

    log.info("finished.")
