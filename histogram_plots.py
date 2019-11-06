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




def histogram(path, data_all,table_path, pdf_name , data_name, xpeak_key, table_key=''):
    pdf_pages = PdfPages(path+pdf_name)
    #print data_path
    #data_all = np.load(data_path)
    tab1 = asc.read(table_path)
    if table_key == 'ID':
        ID = tab1['ID']
        
        calculated_com = tab1['com_O_log']
        calculated_Te = tab1 ['Temperature']
        calculated_single = tab1['O_s_ion']
        calculated_double = tab1['O_d_ion']
        detect = tab1['Detection']
        detection = np.where((detect==1))[0]
        ID_detect = ID[detection]

    '''statistical_all = np.load(statistical_peaks)
    xpeaks = statistical_all[xpeak_key]
    print xpeaks.shape
    print data_hist.shape'''

    for bb in range(len(xpeak_key)):
        for ii in range(len(data_all)):
            print bb, ii
            hist_name = data_name[bb][ii]
            xpeak_name = xpeak_key[ii]
            print hist_name, xpeak_name
            #data_hist = data_all[hist_name]
            data_xpeak = data_all[xpeak_name]
            
            #This is defining a quick fix for passing in several wanted histograms 
            if hist_name == 'Te_propdist': calculated_value = calculated_Te[detection]
            if hist_name == 'O_d_ion_pdf': calculated_value = calculated_double[detection]
            if hist_name == 'O_s_ion_pdf': calculated_value = calculated_single[detection]
            if hist_name == 'com_O_log_pdf': calculated_value = calculated_com[detection]
            
            print(calculated_value)
            
            
            

        nrows = len(calculated_value)/2
        print nrows
        rows =np.arange(nrows)
        print 'a', rows
        ncols = 2
    
    
        for aa in range(xpeaks.shape[0]):
            row = rows[aa/2]
            col = aa% ncols
            #print row, col
            if aa % (nrows*ncols) ==0:
                fig, ax_arr = plt.subplots(nrows = nrows, ncols= ncols, squeeze =False)

            ax = ax_arr[row,col]
            min_val = np.nanmin(data_hist[aa])
            max_val = np.nanmax(data_hist[aa])
            #print min_val, max_val
        
            bin_arr = np.linspace(min_val, max_val)
            #print bin_arr
        
            title ='Bin: ',ID_detect[aa]
            ax.hist(data_hist[aa], bins =bin_arr)
            ax.axvline(x = xpeaks[aa],color = 'r', label = 'compute_one_sig_xpeak', linewidth = 0.5)
            ax.axvline(x = calculated_value[aa],color = 'm', label = 'stacked metallicities', linewidth =0.5)
            ax.set_xlim(min_val, max_val)
            ax.annotate(title, [0.95,0.5], xycoords = 'axes fraction',va = 'center', ha = 'right', fontsize = 6)
            plt.subplots_adjust(left= 0.07 , bottom= 0.10 , right= 0.97, top= 0.97, wspace = 0.15, hspace =0.15)
            if aa%(nrows*ncols) == 0:ax.legend(fontsize = 'xx-small')
            if row != 3:
                ax.set_xticklabels([])
            if row == 3: plt.xlabel("Metallicity") 
            fig.savefig(pdf_pages, format ='pdf')
        
    pdf_pages.close()





def run_histogram_TM(fitspath, TM_file, dict_list, data_name):     

    #dict_list = [Te_pdf_dict,Te_xpeak_dict, metallicity_pdf_dict, metallicity_xpeak_dict]
    #data_name = [['Te_propdist'], ['Te_xpeak'],['O_d_ion_pdf', 'O_s_ion_pdf', 'com_O_log_pdf'],['O_s_ion_xpeak','O_d_ion_xpeak', 'com_O_log_xpeak']]


    #['O_d_ion_pdf', 'O_s_ion_pdf', 'com_O_log_pdf'],
    #['O_d_ion_pdf', 'O_s_ion_pdf', 'com_O_log_pdf']]



    #dictxpeak_key = ['Te_xpeak','O_d_ion_pdf', 'O_s_ion_pdf', 'com_O_log_pdf']




    xpeak_key = ['Te_xpeak','O_s_ion_xpeak','O_d_ion_xpeak', 'com_O_log_xpeak']

    histo_dict = {} #will have all the data and xpeaks for all histograms wanted 
    for bb in range(len(dict_list)):
        dic0 = np.load(dict_list[bb])   ###dictionary = np.load(path of the dictionary)
        #print dic0.keys()
        for aa in range(len(data_name[bb])):
            print data_name[bb]
            print data_name[bb][aa]
            key_name = data_name[bb][aa]     ###''name of the key
            #print key_name
            dic1 = {key_name: dic0[key_name]}             ###{key: dictionary[''name of key]}
            histo_dict.update(dic1)         ###updates new dictionary 
    
    print histo_dict

    
    pdf_name = 'Te_M_histogram_plots.pdf'
    histogram(fitspath, histo_dict, TM_file, pdf_name, data_name, xpeak_key, table_key = 'ID')



    ''''
    for ii in range(len(hist_wanted)):
        #data_name = hist_wanted[ii]
        #print data_name
        if data_name == 'Te_propdist':
            print 'Te'
            data_path = Te_pdf_dict
            print data_path
            statistical_peaks = Te_xpeak_dict
            
        else:
            data_path = metallicity_pdf_dict
            statistical_peaks = metallicity_xpeak_dict
        pdf_name = data_name+'histogram_plots.pdf'
        print pdf_name
        histogram(fitspath, data_path, TM_file,pdf_name, statistical_peaks, table_key = 'ID', data_name = data_name,xpeak_key= xpeak_key[ii])

        histogram(fitspath, histo_dict, TM_file, pdf_name, table_key = 'ID', data_name, xpeak_key)'''



###generalize code further so that the run function passes in a list of stuff
