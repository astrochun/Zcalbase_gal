######Creates histograms of data inputted in dictionaries#####
######Currently set up to plot tempature and metallicity error propagation#####
'''

##############Functions#######################
# histogram(path, data_all, table_path, pdf_name, xpeak_key, table_key='')
Input variables: 
path            -> location of where the outputted pdf_file will be saved
data_all        -> big dictionary of all the error propagation data that will be put into 
                   histogram plot; created in run_histogram_TM
table_path      -> location of the temperature_metallicity table outputted by the R_temp_cal 
                   functions; can also be the combine_flux_table created by 
                   zoom_and_gauss_general
pdf_name        -> name of the outputted pdf file
xpeak_key       -> list of names of the median values outputted from compute_one_sig as xpeak
                   used for naming and referencing
table_key       -> name of one of the columns of the table inputted by table_path
                   used to call the binned data



# run_histogram_TM(path,table_path, dict_list,xpeak_key): 
Input variables: 
path            -> name of where you are working and location of where the 
                   outputted pdf_file will be saved
table_path      -> location of the temperature_metallicity table outputted by the R_temp_cal 
                   functions; can also be the combine_flux_table created by 
                   zoom_and_gauss_general
dict_list       -> list of dictionaries whose data we want to plot in a histogram
xpeak_key       -> list of names of the median values outputted from compute_one_sig as xpeak
                   used for naming and referencing
#Example :
    #         dict_list = [Te_pdf_dict,Te_xpeak_dict,
    #                     metallicity_pdf_dict, metallicity_xpeak_dict]
    #         xpeak_key = ['Te_xpeak','O_s_ion_xpeak',
    #                     'O_d_ion_xpeak', 'com_O_log_xpeak']


Calling order: call run_histogram to combine all dictionaries into one large dictionary that gets passed into histogram and from there is plotted 

'''


import numpy as np
import matplotlib.pyplot as plt
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
from collections import OrderedDict




def histogram(path, data_all,table_path, pdf_name , xpeak_key, table_key=''):
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
        calculated_logs = tab1['log_O_s']
        calculated_logd = tab1['log_O_d']
        detect = tab1['Detection']
        detection = np.where((detect==1))[0]
        ID_detect = ID[detection]

    '''statistical_all = np.load(statistical_peaks)
    xpeaks = statistical_all[xpeak_key]
    print xpeaks.shape
    print data_hist.shape'''

    #print data_name.shape()

    histo_keys = data_all.keys()
    print(histo_keys)
    if type(histo_keys) != list:
        histo_keys = list(histo_keys)

    pdf_list = [histo_keys[xx] for xx in range(len(histo_keys)) if ('pdf' in histo_keys[xx])]
    print(pdf_list)
    xpeak_list = [str0.replace('pdf','xpeak') for str0 in pdf_list]
    print(xpeak_list)






    for ii in range(len(pdf_list)):
            
        hist_name = pdf_list[ii]
        xpeak_name = xpeak_list[ii]
            
        data_hist = data_all[hist_name]
        data_xpeak = data_all[xpeak_name]

        print(xpeak_name)
            
        #This is defining a quick fix for passing in several wanted histograms 
        if hist_name == 'Te_pdf': calculated_value = calculated_Te[detection]
        if hist_name == 'O_d_ion_pdf': calculated_value = calculated_double[detection]
        if hist_name == 'O_s_ion_pdf': calculated_value = calculated_single[detection]
        if hist_name == 'com_O_log_pdf': calculated_value = calculated_com[detection]
        if hist_name == 'O_d_ion_log_pdf': calculated_value = calculated_logd[detection]
        if hist_name == 'O_s_ion_log_pdf': calculated_value = calculated_logs[detection]
        #print('Should all be related:', calculated_value, hist_name, xpeak_name)
            

        print('xpeak: ', data_xpeak)
        print('stacked: ', calculated_value)
            

       
            
    
        if len(calculated_value) % 2 == 0:
            nrows = len(calculated_value)//2
        else:
            nrows = len(calculated_value)//2 + 1

        #print nrows
        rows =np.arange(nrows)
        #print 'a', rows
        ncols = 2
    
    
        for aa in range(len(calculated_value)):
            row = rows[aa//2]
            col = aa% ncols
            #print row, col
            if aa % (nrows*ncols) ==0:
                fig, ax_arr = plt.subplots(nrows = nrows, ncols= ncols, sharex=True, squeeze =False)

            ax = ax_arr[row,col]
            non_inf = np.where(np.isfinite(data_hist[aa]) == True)[0]
            min_val = np.nanmin(data_hist[aa][non_inf])
            max_val = np.nanmax(data_hist[aa][non_inf])
            #print min_val, max_val
        
            bin_arr = np.linspace(min_val, max_val)
            #print bin_arr
        
            title ='Bin: ',ID_detect[aa]
            ax.hist(data_hist[aa][non_inf], bins =bin_arr)
            ax.axvline(x = data_xpeak[aa],color = 'r', label = 'compute_one_sig_xpeak', linewidth = 0.5)
            ax.axvline(x = calculated_value[aa],color = 'm', label = 'stacked metallicities', linewidth =0.5)
            ax.set_xlim(min_val, max_val)
            ax.annotate(title, [0.95,0.5], xycoords = 'axes fraction',va = 'center', ha = 'right', fontsize = 6)
            plt.subplots_adjust(left= 0.07 , bottom= 0.10 , right= 0.97, top= 0.97, wspace = 0.15, hspace =0.15)
            if aa%(nrows*ncols) == 0:
                ax.legend(title = hist_name, fontsize = 3) #fontsize = 'xx-small')
            #if row == 3: plt.xlabel(pdf_list[ii])
        fig.savefig(pdf_pages, format ='pdf')
        
    pdf_pages.close()





def run_histogram(path, table_path, dict_list, xpeak_key):    #,data_name):
    if path[-1] != "/": path +="/"
    histo_dict = OrderedDict()  #will have all the data and xpeaks for all histograms wanted 
    for bb in range(len(dict_list)):
        dic0 = np.load(dict_list[bb])   ###dictionary = np.load(path of the dictionary)
        dic0_keys = dic0.keys()
        print(dic0_keys)
        histo_dict.update(dic0)
        '''for key in dict0_keys:
            dic1 = {key: dic0[key]}             ###{key: dictionary[''name of key]}
            histo_dict.update(dic1)'''       ###updates new dictionary 
    
    print(histo_dict)

    
    pdf_name = 'Te_M_histogram_plots.pdf'
    histogram(fitspath, histo_dict, table_path, pdf_name, xpeak_key, table_key = 'ID')



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



