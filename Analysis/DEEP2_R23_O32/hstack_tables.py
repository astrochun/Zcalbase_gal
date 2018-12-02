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
import glob


#This code stacks ascii flux tables generated by the zoom_and_gauss_general code into a large ascii file


def h_stack(fitspath, table_files, asc_intro,new_name):
    table_files.sort()
    print 'import'
    for ii in range(len(table_files)):
        asc_tab= asc.read(table_files[ii])
        if ii == 0: 
            hstacking = hstack([asc_intro,asc_tab])
        else:
            hstacking = hstack([hstacking,asc_tab])

        print type(hstacking)
    #Grid Fits Files
    #hstacking.write(new_name, format='fits')
    #Grid Ascii Files
    asc.write(hstacking,new_name, overwrite = True)













