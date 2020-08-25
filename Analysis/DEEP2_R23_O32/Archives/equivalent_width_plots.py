
"""
PLOTTING RESULTS FROM ANALYSIS
Some were written before MSC and are now in MSC

Functions:
-Equivalent Width for R23
-Equivalent Width for O32
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import zoom_and_gauss_general
from Metallicity_Stack_Commons.Metallicity_Stack_Commons import lambda0, line_type, line_name

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'


def ew_plot_r23(fitspath, asc_table):
    """
    Purpose
    Plots equivalent width verse R23 and saves to a pdf file.
    """
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
            ax.set_title('EW vs. R23   ' + str(np.int(lambda0[oo])))
            fig.set_size_inches(8, 8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()
    pdf_pages.close()


def ew_plot_o32(fitspath, asc_table):
    """
    Purpose
    Plots equivalent width verse O32 and saves to a pdf file.
    """
    pdf_pages = PdfPages(join(fitspath, 'equivalent_vs_O32_plots.pdf'))
    asc_tab = asc.read(asc_table)
    O32 = asc_tab['logO32_avg'].data
    for oo in range(len(lambda0)):
        if line_type[oo] == 'Balmer':
            equ_w = asc_tab['EW_'+str(np.int(lambda0[oo]))+'_abs'].data
            fig, ax = plt.subplots()
            ax.scatter(equ_w, O32, marker='.')
            ax.set_xlabel('Equivalent Width')
            ax.set_ylabel('O32')
            ax.set_title('EW vs. O32   ' + str(np.int(lambda0[oo])))
            fig.set_size_inches(8, 8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()
    pdf_pages.close()
