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
table_key       -> name of one of the columns of the table inputted by table_path
                   used to call the binned data

Output: pdf file with all the histograms plotted  

# run_histogram_TM(path,table_path, dict_list,xpeak_key): 
Input variables: 
path            -> name of where you are working and location of where the 
                   outputted pdf_file will be saved
table_path      -> location of the temperature_metallicity table outputted by the R_temp_cal 
                   functions; can also be the combine_flux_table created by 
                   zoom_and_gauss_general
dict_list       -> list of dictionaries whose data we want to plot in a histogram

#Example :
    #         dict_list = [Te_pdf_dict,Te_xpeak_dict,
    #                     metallicity_pdf_dict, metallicity_xpeak_dict]
    #         xpeak_key = ['Te_xpeak','O_s_ion_xpeak',
    #                     'O_d_ion_xpeak', 'com_O_log_xpeak']


Calling order: call run_histogram to combine all dictionaries into one large dictionary that gets passed into histogram and from there is plotted 



#run_histogram_FR(path, table_path, dict_list, verification_table, sharex = False)

#Notes: This function takes the flux_pdf distributions computed in error_prop and calculates the flux ratios. Then it calls compute_onesig_pdf from chun_codes for the flux ratios. Finally it create the compute dictionary of ratios, xpeaks, and low/high limits and passes it all into histogram
Input Variables 
path            -> name of where you are working and location of where the 
                   outputted pdf_file will be saved
table_path      -> location of the temperature_metallicity table outputted by the R_temp_cal 
                   functions; can also be the combine_flux_table created by 
                   zoom_and_gauss_general
dict_list       -> list of dictionaries whose data we want to plot in a histogram
verification_table -> table created independently of this code that tells us whether or not each bin is a detection

Outputs: pdf file of all the histograms
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
from datetime import date

from chun_codes import compute_onesig_pdf


def histogram(path, data_all,table_path, pdf_name,  verification_table, table_key='', sharex = False):
    
    pdf_pages = PdfPages(path+pdf_name)
    #print data_path
    #data_all = np.load(data_path)
    tab1 = asc.read(table_path)
    verify = asc.read(verification_table)
    detect = verify['Detection']
    ID = verify['ID']
    detection = np.where((detect==1))[0]
    ID_detect = ID[detection]
    print(ID_detect)

    #For plotting Temperature and Metallicity Histograms 
    if table_key == 'Temperature':
        #ID = tab1['ID']
        
        calculated_com = tab1['com_O_log']
        calculated_Te = tab1 ['Temperature']
        calculated_single = tab1['O_s_ion']
        calculated_double = tab1['O_d_ion']
        calculated_logs = tab1['log_O_s']
        calculated_logd = tab1['log_O_d']
        

        
       
            
    #For plotting Flux Ratio Histograms ['OIII_4363','OII_3727','OIII_4958','HGAMMA','OIII_5007','HBETA','HDELTA']
    if table_key == 'ID':
        O5007 = tab1['OIII_5007_Flux_Observed']
        O4363 = tab1 ['OIII_4363_Flux_Observed']
        O4959 = tab1['OIII_4958_Flux_Observed']
        HBETA = tab1['HBETA_Flux_Observed']
        O3727 = tab1['OII_3727_Flux_Observed']

        O3727_HBETA = O3727[detection]/HBETA[detection]
        O5007_HBETA = O5007[detection]/HBETA[detection]
        #O4959_HBETA = O4959[detection]/HBETA[detection]
        O4363_O5007  = O4363[detection]/O5007[detection]
        O5007_O4958  = O5007[detection]/O4959[detection]
        #O4363_O4959  = O4363[detection]/O4959[detection]
    
        O5007_O3727  = O5007[detection]/O3727[detection]
        #O4959_O3727  = O4959[detection]/O3727[detection]

        R23_combine = (O3727_HBETA + O5007_HBETA*(1+1/3.1))/HBETA[detection]

        #print('O3727_HBETA', O3727_HBETA, 'O5007_HBETA', O5007_HBETA, 'O4959_HBETA', O4959_HBETA, 'O4363_O5007', O4363_O5007, 'O4363_O4959', O4363_O4959, 'O5007_O3727', O5007_O3727, 'O4959_O3727', O4959_O3727)
        

    histo_keys = data_all.keys()
    #print(histo_keys)
    if type(histo_keys) != list:
        histo_keys = list(histo_keys)
            
    xpeak_list = [histo_keys[xx] for xx in range(len(histo_keys)) if ('xpeak' in histo_keys[xx])]
    pdf_list = [str0.replace('_xpeak','' ) for str0 in xpeak_list]
    print(xpeak_list)
    error_list = [str0.replace('xpeak','lowhigh_error') for str0 in pdf_list]
    print(error_list)
    print(len(pdf_list), pdf_list)





    for ii in range(len(pdf_list)):
            
        hist_name = pdf_list[ii]
        print(hist_name)
        xpeak_name = xpeak_list[ii]
        error_name = error_list[ii]
            
        data_hist = data_all[hist_name]
        data_xpeak = data_all[xpeak_name]
        data_errors = data_all[error_name]

        data_lowerror = data_errors[:,0]
        data_higherror= data_errors[:,1]
        #print(data_errors)
        #print('Low: ',data_lowerror)
        #print('High: ',data_higherror)
        
            
        #This is defining a quick fix for passing in several wanted histograms 
        if hist_name == 'Te_pdf': calculated_value = calculated_Te[detection]
        if hist_name == 'O_d_ion_pdf': calculated_value = calculated_double[detection]
        if hist_name == 'O_s_ion_pdf': calculated_value = calculated_single[detection]
        if hist_name == 'com_O_log_pdf': calculated_value = calculated_com[detection]
        if hist_name == 'O_d_ion_log_pdf': calculated_value = calculated_logd[detection]
        if hist_name == 'O_s_ion_log_pdf': calculated_value = calculated_logs[detection]

        
        if hist_name == '3727/HBETA': calculated_value = O3727_HBETA
        if hist_name == '5007/HBETA': calculated_value = O5007_HBETA
        #if hist_name == '4959/HBETA': calculated_value = O4959_HBETA
        if hist_name == '4363/5007': calculated_value = O4363_O5007
        #if hist_name == '4363/4959': calculated_value =  O4363_O4959
        if hist_name == '5007/3727': calculated_value = O5007_O3727
        #if hist_name == '4959/3727' : calculated_value = O4959_O3727
        if hist_name == '5007/4959': calculated_value = O5007_O4958
        if hist_name == 'R23': calculated_value = R23_combine
 
    
        if len(calculated_value) % 2 == 0:
            nrows = len(calculated_value)//2
        else:
            nrows = len(calculated_value)//2 + 1

        print('nrows', nrows)
        rows =np.arange(nrows)
        print 'rows', rows
        ncols = 2
    
    
        for aa in range(len(calculated_value)):
            print(aa)
            row = rows[aa//2]
            col = aa% ncols
            #print row, col
            if aa % (nrows*ncols) ==0:
                fig, ax_arr = plt.subplots(nrows = nrows, ncols= ncols,  sharex=sharex,squeeze =False)  

            ax = ax_arr[row,col]
            print(row,col)

            non_inf = np.where(np.isfinite(data_hist[aa]) == True)[0]
            min_val = np.nanmin(data_hist[aa][non_inf])
            max_val = np.nanmax(data_hist[aa][non_inf])
            
            #print(data_xpeak[aa], data_lowerror[aa], data_higherror[aa])
            lowerr= data_xpeak[aa]-data_lowerror[aa]
            higherr =  data_xpeak[aa]+data_higherror[aa]
            
            
            bin_arr = np.linspace(min_val, max_val)
            
            title ='Bin: ',ID_detect[aa]
            ax.hist(data_hist[aa][non_inf], bins =bin_arr)
            ax.axvline(x = data_xpeak[aa],color = 'r', label = 'compute_one_sig_xpeak', linewidth = 0.5)
            ax.axvline(x = calculated_value[aa],color = 'm', label = 'stacked metallicities', linewidth =0.5)
            #ax.axvline(x = lowerr,color = 'g', label = 'Lower Limit', linewidth =0.5)
            #ax.axvline(x = higherr,color = 'k', label = 'Upper Limit', linewidth =0.5)
            


            #ax.set_xlim(min_val, max_val)
            ax.annotate(title, [0.95,0.5], xycoords = 'axes fraction',va = 'center', ha = 'right', fontsize = 6)
            if sharex == True: plt.subplots_adjust(left= 0.07 , bottom= 0.10 , right= 0.97, top= 0.97, wspace = 0.15, hspace =0.15)
            if sharex == False:
                plt.subplots_adjust(left= 0.07 , bottom= 0.05 , right= 0.97, top= 0.97, wspace = 0.15, hspace =0.2)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(7) 
            if aa%(nrows*ncols) == 0:
                ax.legend(title = hist_name, fontsize = 3) 
                #if row == 3: plt.xlabel(pdf_list[ii])
        fig.savefig(pdf_pages, format ='pdf')
        
    pdf_pages.close()





#######################Calling Functions###########################


###Temperature and Metallicity Dictionary List for Reagen###
#dict_list = [Te_pdf_dict,Te_xpeak_dict, metallicity_pdf_dict, metallicity_xpeak_dict, Te_error_dict, metallicity_error_dict]

###Flux Dictionary List for Reagen###
#dict_list = [flux_pdf_dict, flux_errors, flux_lowhigh]

def run_histogram_TM(path, table_path, dict_list, sharex = False):    #,data_name):
    if path[-1] != "/": path +="/"
    histo_dict = OrderedDict()  #will have all the data and xpeaks for all histograms wanted 
    for bb in range(len(dict_list)):
        dic0 = np.load(dict_list[bb])   ###dictionary = np.load(path of the dictionary)
        #dic0_keys = dic0.keys()
        #print(dic0_keys)
        histo_dict.update(dic0)
        '''for key in dict0_keys:
            dic1 = {key: dic0[key]}             ###{key: dictionary[''name of key]}
            histo_dict.update(dic1)'''       ###updates new dictionary 
    
    #print(histo_dict)

    today = date.today()
    pdf_name = 'Te_M_histogram_plots_' + "%02i%02i" % (today.month, today.day)+'.pdf'
    histogram(path, histo_dict, table_path, pdf_name, verification_table, table_key = 'Temperature', sharex = sharex)




###Flux Dictionary List for Reagen###
#dict_list = [flux_pdf_dict, flux_errors, flux_lowhigh]



#This call will create histograms for the flux ratios
#We need this second call because the flux dictionaries have to be changed into ratios

#This call calculates the flux ratios from the combine_flux_ascii table and random_pdf (of each line), gets the one_sig_pdf measurements for the ratios, and saves all in dictionary before 
def run_histogram_FR(path, table_path, dict_list, verification_table, sharex = False):    #,data_name):
    if path[-1] != "/": path +="/"

    tab1 = asc.read(table_path)
    O5007 = tab1['OIII_5007_Flux_Observed']
    O4363 = tab1 ['OIII_4363_Flux_Observed']
    O4959 = tab1['OIII_4958_Flux_Observed']
    HBETA = tab1['HBETA_Flux_Observed']
    O3727 = tab1['OII_3727_Flux_Observed']

    O3727_HBETA = O3727/HBETA
    O5007_HBETA = O5007/HBETA
    O4959_HBETA = O4959/HBETA
    O4363_O5007  = O4363/O5007
    O4363_O4959  = O4363/O4959
    
    O5007_O3727  = O5007/O3727
    O4959_O3727  = O4959/O3727

    O5007_O4958 = O5007/O4959
    R23_combine = (O3727_HBETA+O5007_HBETA*(1+1/3.1))/HBETA
    
    flux_dict = np.load(dict_list[0])
    
    #Two_Beta = OII/HBETA
    #Three_Beta = (OIII5007+OIII4959)/HBETA
    ####Making the flux_ratio dictionary########
    c3727_HBETA = flux_dict['OII_3727'] / flux_dict['HBETA'] 
    c5007_HBETA = flux_dict['OIII_5007']/ flux_dict['HBETA'] 
    #c4959_HBETA = flux_dict['OIII_4958']/ flux_dict['HBETA']
    c4363_c5007 = flux_dict['OIII_4363']/ flux_dict['OIII_5007']
    c5007_c4958 = flux_dict['OIII_5007']/ flux_dict['OIII_4958']
    #c4363_c4959 = flux_dict['OIII_4363']/ flux_dict['OIII_4958']
    R23 = (c3727_HBETA + c5007_HBETA*(1+1/3.1)) /flux_dict['HBETA']
    
    c5007_c3727  = flux_dict['OIII_5007']/ flux_dict['OII_3727']
    #c4959_c3727  = flux_dict['OIII_4958']/ flux_dict['OII_3727']
        

    pdf_list =  ['3727/HBETA', '5007/HBETA', '4363/5007', '5007/4959',  '5007/3727', 'R23']
    data_list = [c3727_HBETA,   c5007_HBETA,  c4363_c5007, c5007_c4958,  c5007_c3727, R23]
    ratio_list =[O3727_HBETA,   O5007_HBETA , O4363_O5007, O5007_O4958 , O5007_O3727, R23_combine]

    flux_ratio_dict = {}
    fr_xpeak_dict = {}
    fr_error_dict = {}
    
    for ii in range(len(pdf_list)):
        flux_ratio_dict[pdf_list[ii]] = data_list[ii]

        error, xpeak = compute_onesig_pdf(data_list[ii], ratio_list[ii], usepeak=True, silent=True, verbose = True)
        #print('EX:', error, xpeak)
        fr_xpeak_dict[pdf_list[ii]+'_xpeak'] = xpeak
        fr_error_dict[pdf_list[ii]+'_errors'] = error

    new_dict_list = [flux_ratio_dict, fr_xpeak_dict, fr_error_dict]
    histo_dict = OrderedDict()  #will have all the data and xpeaks for all histograms wanted 
    for bb in range(len(new_dict_list)):
        dic0 = new_dict_list[bb]   ###dictionary = np.load(path of the dictionary)
        
        histo_dict.update(dic0)
        '''for key in dict0_keys:
            dic1 = {key: dic0[key]}             ###{key: dictionary[''name of key]}
            histo_dict.update(dic1)'''       ###updates new dictionary 
    
    #print(histo_dict)

    
    pdf_name = 'Flux_Ratio_histogram_plots.pdf'
    histogram(path, histo_dict, table_path, pdf_name, verification_table, table_key = 'ID', sharex = sharex)


