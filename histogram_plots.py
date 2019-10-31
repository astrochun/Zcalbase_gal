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
from mpl_toolkits.mplot3d import Axes3D
import sys




def histogram(path, data_path,table_path, pdf_name , statistical_peaks, table_key='', data_name = ''): 
    data_all = np.load(data_path)
    data_hist = data_all[data_name]

    statistical_all = np.load(statistical_peaks)
    xpeaks = statistical_all[data_name]
    print xpeaks.shape[0]

    tab1 = asc.read(table_path)
    if table_key == 'ID':
        ID = tab1['ID']
        
        calculated_metal = tab1['com_O_log']
        detect = tab1['Detection']
        detection = np.where((detect==1))[0]
        com_metal = calculated_metal[detection]
        ID_detect = ID[detection]
    pdf_pages = PdfPages(path+pdf_name)
    
    
    
    for aa in range(xpeaks.shape[0]):
        fig, ax = plt.subplots()
        min_val = np.nanmin(data_hist[aa])
        max_val = np.nanmax(data_hist[aa])
        print min_val, max_val
        
        bin_arr = np.linspace(min_val, max_val)
        #print bin_arr
        
        
        ax.hist(data_hist[aa], bins =bin_arr)
        ax.axvline(x = xpeaks[aa],color = 'r', label = 'compute_one_sig_xpeak')
        ax.axvline(x = com_metal[aa],color = 'm', label = 'stacked metallicities')
        title ='Bin: ',ID_detect[aa]     #,'Detection: ',detect[aa]
        #print title
        ax.set_title(title)
        ax.legend()
        fig.savefig(pdf_pages, format ='pdf')
        fig.clear()
    pdf_pages.close()
