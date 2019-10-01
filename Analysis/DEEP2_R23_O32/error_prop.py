

import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf


#Calling Error Propogation
#error_prop_chuncodes(fitspath, dataset, working_wave, flux_g_array, RMS_array)

def error_prop_chuncodes(fitspath, dataset, working_wave, values,RMS):
    #x_arr = np.random.random_integers(10,10)
    #all emission lines
    #p_page= fitspath+dataset+'/'+str(np.int(working_wave))+'error_histogram.pdf'
    #print p_page
    #pdf_pages_err = PdfPages(p_page)
    #print 'values:', values
    #print 'RMS:', RMS
    fluxg_pdf = random_pdf(values,RMS, seed_i =1, n_iter=1000, silent = False)
    print 'err_function:', fluxg_pdf
    print 'ran errors'
    #fig_err, ax = plt.subplots()
    #plt.hist(fluxg_pdf, 50)

    #fig_err.set_size_inches(8,8)
    #fig_err.savefig(pdf_pages_err, format='pdf')
    #pdf_pages_err.close()
    
    return fluxg_pdf
