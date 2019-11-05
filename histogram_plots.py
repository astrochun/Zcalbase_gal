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
    print xpeaks.shape
    print data_hist.shape

    tab1 = asc.read(table_path)
    if table_key == 'ID':
        ID = tab1['ID']
        
        calculated_metal = tab1['com_O_log']
        detect = tab1['Detection']
        detection = np.where((detect==1))[0]
        com_metal = calculated_metal[detection]
        ID_detect = ID[detection]


    pdf_pages = PdfPages(path+pdf_name)

    nrows = len(com_metal)/2
    print nrows
    rows =np.arange(nrows)
    print 'a', rows
    ncols = 2
    
    
    for aa in range(xpeaks.shape[0]):
        row = rows[aa/2]
        col = aa% ncols
        print row, col
        if aa % (nrows*ncols) ==0:
            fig, ax_arr = plt.subplots(nrows = nrows, ncols= ncols, squeeze =False)

        ax = ax_arr[row,col]
        min_val = np.nanmin(data_hist[aa])
        max_val = np.nanmax(data_hist[aa])
        print min_val, max_val
        
        bin_arr = np.linspace(min_val, max_val)
        #print bin_arr
        
        title ='Bin: ',ID_detect[aa]
        ax.hist(data_hist[aa], bins =bin_arr)
        ax.axvline(x = xpeaks[aa],color = 'r', label = 'compute_one_sig_xpeak', linewidth = 0.5)
        ax.axvline(x = com_metal[aa],color = 'm', label = 'stacked metallicities', linewidth =0.5)
        ax.set_xlim(7.75,9.15)
        ax.annotate(title, [0.95,0.5], xycoords = 'axes fraction',va = 'center', ha = 'right', fontsize = 6)
        plt.subplots_adjust(left= 0.05 , bottom= 0.10 , right= 0.95, top= 0.90, wspace = 0.15, hspace =0.15)
        if aa%(nrows*ncols) == 0:ax.legend(fontsize = 'xx-small')
        if row != 3:
            ax.set_xticklabels([])
        if row == 3: plt.xlabel("Metallicity") 
    fig.savefig(pdf_pages, format ='pdf')
        
    pdf_pages.close()





def run_histogram_TM(fitspath,Te_pdf_dict,Te_xpeak_dict, metallicity_pdf_dict, metallicity_xpeak_dict, TM_file):

    hist_wanted = ['Te_propdist', 'com_O_log,', 'O_s_ion', 'O_d_ion']

    for ii in range(len(hist_wanted)):
        data_name = hist_wanted[ii]
        print data_name
        if data_name == 'Te':
            data_path = Te_pdf_dict
            statistical_peaks = Te_xpeak_dict
            
        else:
            data_path = metallicity_pdf_dict
            statistical_peaks = metallicity_xpeak_dict
        pdf_name = data_name+'histogram_plots.pdf'
        histogram(fitspath, data_path, TM_file,pdf_name, statistical_peaks, table_key = 'ID', data_name = data_name)           
