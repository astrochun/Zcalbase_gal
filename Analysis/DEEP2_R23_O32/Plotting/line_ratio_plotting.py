
###CREATES LINE PLOTS FOR CHECKING RESULTS
###CALLED EVERYTIME GENERAL FUNCTIONS ARE CALLED

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

import general


def Plotting_Data1(fitspath, dataset, combine_flux_ascii, asc_table1):
    
    line_plot = fitspath +dataset+'_line_ratio_plots.pdf'
    
    #combine_flux_ascii = fitspath+dataset#+'_combined_flux_table.tbl'
    print "### combine_flux_ascii : "+combine_flux_ascii 
    fitted_data = asc.read(combine_flux_ascii)

    print "### asc_table1 : "+asc_table1
    raw_data = asc.read(asc_table1)
    
    OII = fitted_data['OII_3727_Flux_Observed']
    OIII4959 = fitted_data['OIII_4958_Flux_Observed']
    OIII5007 = fitted_data['OIII_5007_Flux_Observed']
    H_BETA = fitted_data['HBETA_Flux_Observed']
    binnum = fitted_data['N_Galaxies']
    ID = fitted_data['ID']
    print 'binnum:', binnum, len(binnum)
    pdf_pages = PdfPages(line_plot)
    nrows = 4
    ncols = 4

    R23_composite = np.zeros(binnum.shape[0])
    O32_composite = np.zeros(binnum.shape[0]) 
    for ii in range(len(binnum)):
        R23_comp = np.log10((OII[ii] + (1.33*OIII5007[ii]))/H_BETA[ii])
        O32_comp = np.log10((1.33*OIII5007[ii])/OII[ii])
        print R23_comp, O32_comp
        R23_composite[ii]= R23_comp
        O32_composite[ii]= O32_comp
    
    R23_raw = raw_data['xBar']
    O32_raw = raw_data['yBar']
    binnum_raw = raw_data['area']

    if dataset != 'Grid':
        for rr in range(len(binnum)):
            if binnum[rr] == binnum_raw[rr]:
                print 'equal',binnum[rr], binnum_raw[rr]

    print 'R23_raw: as calculated by the grid or voronoi code', R23_raw, 'O32_raw: as calculated by the grid or voronoi code', O32_raw
    print 'R23_composite: as calculated from observations', R23_composite, 'O32_composite: as calculated from observations', O32_composite 
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(R23_raw,R23_composite, marker= 'o', facecolor= 'none', edgecolor ='b',label= 'R23 Ratio: Vornoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset+' Raw vs. Composite for R23')
    for rr in range(len(ID)):
        ax_arr.annotate(ID[rr], (R23_raw[rr], R23_composite[rr]))
    ax_arr.set_xlabel(r'Raw log($R_{23}$)')
    ax_arr.set_ylabel(r'Composite log($R_{23}$)')
    ax_arr.plot([0.0,1.3], [0.0,1.3], 'k-')
    
    fig.savefig(pdf_pages, format='pdf')

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(O32_raw,O32_composite, marker= 'o', facecolor= 'none', edgecolor ='b', label= 'O32 Ratio: Vornoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset+'Raw vs. Composite for O32')
    for oo in range(len(ID)):
        ax_arr.annotate(ID[oo], (O32_raw[oo], O32_composite[oo]))
    ax_arr.set_xlabel(r'Raw log($O_{32}$)')
    ax_arr.set_ylabel(r'Composite log($O_{32}$)')

    ax_arr.plot([-1,1.2], [-1,1.2], 'k-')
    fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()

    fig.clear()
