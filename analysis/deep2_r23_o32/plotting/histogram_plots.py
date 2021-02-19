# Creates histograms of data inputted in dictionaries
# Currently set up to plot tempature and metallicity error propagation

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
from datetime import date
from os.path import join

from chun_codes import compute_onesig_pdf
from Metallicity_Stack_Commons.column_names import temp_metal_names0, \
    bin_names0


def histogram(path, data_all, input_file, pdf_file,
              valid_file, table_key='', sharex=False):
    """
    Plots histograms for inputted values

    :param path: str. location of where the outputted pdf_file will be saved
    :param data_all: dict. dictionary of all the error propagation data that
                     will be put into histogram plot; created in run_histogram_TM
    :param input_file: str. location of the temperature_metallicity table
                       outputted functions; can also be the
                       combine_flux_table created zoom_and_gauss_general
    :param pdf_file: str. name of the outputted pdf file
    :param table_key: keyword. name of one of the columns of the table inputted
                        by input_file used to call the binned data

    PDF File: path + pdf_file
    No returns
    """
    pp = PdfPages(join(path, pdf_file))
    input_tab = asc.read(input_file)

    # Importing verification table that is always used
    # to get the values for just the detections
    valid_tab = asc.read(valid_file)
    detect = valid_tab[bin_names0[2]]
    ID = valid_tab[bin_names0[0]]
    detection = np.where((detect == 1))[0]
    ID_detect = ID[detection]

    # Importing calculated values from our stacked measures
    # For plotting Temperature and Metallicity Histograms
    if table_key == 'T_e':
        calculated_com = input_tab[temp_metal_names0[1]]
        calculated_Te = input_tab[temp_metal_names0[0]]
        calculated_single = input_tab[temp_metal_names0[4]]
        calculated_double = input_tab[temp_metal_names0[5]]
        calculated_logs = input_tab[temp_metal_names0[2]]
        calculated_logd = input_tab[temp_metal_names0[3]]

    # For plotting Flux Ratio Histograms
    if table_key == 'bin_ID':
        O5007 = input_tab['OIII_5007_Flux_Observed']
        O4363 = input_tab['OIII_4363_Flux_Observed']
        O4959 = input_tab['OIII_4958_Flux_Observed']
        HBETA = input_tab['HBETA_Flux_Observed']
        O3727 = input_tab['OII_3727_Flux_Observed']

        O3727_HBETA = O3727/HBETA
        O5007_HBETA = O5007/HBETA
        O4363_O5007 = O4363/O5007
        O5007_O4958 = O5007/O4959
        O5007_O3727 = O5007/O3727
        R23_combine = (O3727_HBETA + O5007_HBETA * (1 + 1/3.1))

    histo_keys = data_all.keys()
    if type(histo_keys) != list:
        histo_keys = list(histo_keys)

    # Making our lists we will use for plotting using the histogram keys
    xpeak_list = [histo_keys[xx] for xx in
                  range(len(histo_keys)) if ('xpeak' in histo_keys[xx])]
    pdf_list = [str0.replace('_xpeak', '_pdf') for str0 in xpeak_list]
    error_list = [str0.replace('xpeak', 'lowhigh_error')
                  for str0 in xpeak_list]

    # For loop organizes data then starts the histogram plotting
    for ii in range(len(pdf_list)):
        hist_name = pdf_list[ii]
        xpeak_name = xpeak_list[ii]
        error_name = error_list[ii]
            
        data_hist = data_all[hist_name]
        data_xpeak = data_all[xpeak_name]
        data_errors = data_all[error_name]

        data_lowerror = data_errors[:, 0]
        data_higherror = data_errors[:, 1]

        # This is defining a quick fix for passing in several wanted histograms
        if hist_name == 'Te_pdf':
            calculated_value = calculated_Te[detection]
        if hist_name == 'O_d_ion_pdf':
            calculated_value = calculated_double[detection]
        if hist_name == 'O_s_ion_pdf':
            calculated_value = calculated_single[detection]
        if hist_name == 'com_O_log_pdf':
            calculated_value = calculated_com[detection]
        if hist_name == 'O_d_ion_log_pdf':
            calculated_value = calculated_logd[detection]
        if hist_name == 'O_s_ion_log_pdf':
            calculated_value = calculated_logs[detection]

        if hist_name == '3727/HBETA_pdf':
            calculated_value = O3727_HBETA[detection]
        if hist_name == '5007/HBETA_pdf':
            calculated_value = O5007_HBETA[detection]
        if hist_name == '4363/5007_pdf':
            calculated_value = O4363_O5007[detection]
        if hist_name == '5007/3727_pdf':
            calculated_value = O5007_O3727[detection]
        if hist_name == '5007/4959_pdf':
            calculated_value = O5007_O4958[detection]
        if hist_name == 'R23_pdf':
            calculated_value = R23_combine[detection]

        if len(calculated_value) % 2 == 0:
            nrows = len(calculated_value)//2
        else:
            nrows = len(calculated_value)//2 + 1

        print('nrows', nrows)
        rows = np.arange(nrows)
        print('rows', rows)
        ncols = 2

        for aa in range(len(calculated_value)):
            print(aa)
            row = rows[aa//2]
            col = aa % ncols
            if aa % (nrows * ncols) == 0:
                fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols,
                                           sharex=sharex, squeeze=False)

            ax = ax_arr[row, col]
            print(row, col)

            non_inf = np.where(np.isfinite(data_hist[aa]) == True)[0]
            min_val = np.nanmin(data_hist[aa][non_inf])
            max_val = np.nanmax(data_hist[aa][non_inf])

            lowerr = data_xpeak[aa]-data_lowerror[aa]
            higherr = data_xpeak[aa]+data_higherror[aa]

            bin_arr = np.linspace(min_val, max_val)
            
            title = 'Bin: ', ID_detect[aa]
            ax.hist(data_hist[aa][non_inf], bins=bin_arr)
            ax.axvline(x=data_xpeak[aa], color='r',
                       label='compute_one_sig_xpeak', linewidth=0.5)
            ax.axvline(x=calculated_value[aa], color='m',
                       label='stacked metallicities', linewidth=0.5)
            # ax.axvline(x=lowerr, color = 'g',
            #            label='Lower Limit', linewidth=0.5)
            # ax.axvline(x=higherr, color = 'k',
            #            label='Upper Limit', linewidth=0.5)

            # ax.set_xlim(min_val, max_val)
            ax.annotate(title, [0.95, 0.5], xycoords='axes fraction',
                        va='center', ha='right', fontsize=6)
            if sharex:
                plt.subplots_adjust(left=0.07, bottom=0.10, right=0.97,
                                    top=0.97, wspace=0.15, hspace=0.15)
            else:
                plt.subplots_adjust(left=0.07, bottom=0.05, right=0.97,
                                    top=0.97, wspace=0.15, hspace=0.2)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(7) 
            if aa % (nrows * ncols) == 0:
                ax.legend(title=hist_name, fontsize=3)
                # if row == 3: plt.xlabel(pdf_list[ii])
        fig.savefig(pp, format='pdf')
        
    pp.close()


def run_histogram_tm(path, input_file, dict_list,
                     valid_file, sharex=False):
    """
    Call run_histogram to combine all dictionaries into one large dictionary
    that gets passed into histogram and from there is plotted

    :param path: str. name of where you are working and location of where the
                outputted pdf_file will be saved
    :param input_file: str. location of the temperature_metallicity table
                       outputted by the error_prop functions;
                       can also be the combine_flux_table created by
                        zoom_and_gauss_general
    :param dict_list: ls. list of dictionaries whose data we want to plot
                      in a histogram

    Example :
             dict_list = [Te_pdf_dict,Te_xpeak_dict,
                         metallicity_pdf_dict, metallicity_xpeak_dict]
             xpeak_key = ['Te_xpeak','O_s_ion_xpeak',
                         'O_d_ion_xpeak', 'com_O_log_xpeak']

    No Returns
    """
    if path[-1] != "/":
        path += "/"
    # will have all the data and xpeaks for all histograms wanted
    histo_dict = OrderedDict()
    for bb in range(len(dict_list)):
        # dictionary = np.load(path of the dictionary)
        dic0 = np.load(dict_list[bb])
        # Adding small dictionary into big dictionary
        histo_dict.update(dic0)

    today = date.today()
    pdf_file = 'Te_M_histogram_plots_' + "%02i%02i" \
               % (today.month, today.day) + '.pdf'
    histogram(path, histo_dict, input_file, pdf_file, valid_file,
              table_key='T_e', sharex=sharex)


def run_histogram_fr(path, input_file, dict_list, valid_file,
                     sharex=False):
    """
    This function takes the flux_pdf distributions computed in error_prop
    and calculates the flux ratios.
    Then it calls compute_onesig_pdf from chun_codes for the flux ratios.
    Finally it create the compute dictionary of ratios, xpeaks,
    and low/high limits and passes it all into histogram

    :param path: str. name of where you are working and location of where the
                        outputted pdf_file will be saved
    :param input_file: str. location of the flux_ratio table outputted by
                            the R_temp_cal functions; can also be
                            combine_flux_table created zoom_and_gauss_general
    :param dict_list: ls. list of dictionaries whose data we want to plot
                            in a histogram
    :param valid_file: str. table created independently of this code
                        that tells us whether or not each bin is a detection

    No returns
    """
    if path[-1] != "/":
        path += "/"

    valid_tab = asc.read(valid_file)
    detect = valid_tab['Detection']
    detection = np.where((detect == 1))[0]
    
    input_tab = asc.read(input_file)
    tab1 = input_tab[detection]
    O5007 = tab1['OIII_5007_Flux_Observed']
    O4363 = tab1['OIII_4363_Flux_Observed']
    O4959 = tab1['OIII_4958_Flux_Observed']
    HBETA = tab1['HBETA_Flux_Observed']
    O3727 = tab1['OII_3727_Flux_Observed']

    O3727_HBETA = O3727/HBETA
    O5007_HBETA = O5007/HBETA
    O4363_O5007 = O4363/O5007
    
    O5007_O3727 = O5007/O3727

    O5007_O4958 = O5007/O4959
    R23_combine = (O3727_HBETA + O5007_HBETA * (1 + 1/3.1))
    
    flux_dict = np.load(dict_list[0])
    
    # Making the flux_ratio dictionary
    c3727_HBETA = flux_dict['OII_3727_pdf'] / flux_dict['HBETA_pdf'] 
    c5007_HBETA = flux_dict['OIII_5007_pdf'] / flux_dict['HBETA_pdf']
    c4363_c5007 = flux_dict['OIII_4363_pdf'] / flux_dict['OIII_5007_pdf']
    c5007_c4958 = flux_dict['OIII_5007_pdf'] / flux_dict['OIII_4958_pdf']
    R23 = (c3727_HBETA + c5007_HBETA * (1 + 1/3.1))
    c5007_c3727 = flux_dict['OIII_5007_pdf'] / flux_dict['OII_3727_pdf']

    pdf_list = ['3727/HBETA', '5007/HBETA', '4363/5007', '5007/4959',
                '5007/3727', 'R23']
    data_list = [c3727_HBETA, c5007_HBETA, c4363_c5007, c5007_c4958,
                 c5007_c3727, R23]
    ratio_list = [O3727_HBETA, O5007_HBETA, O4363_O5007, O5007_O4958,
                  O5007_O3727, R23_combine]

    flux_ratio_dict = {}
    fr_xpeak_dict = {}
    fr_error_dict = {}
    
    for ii in range(len(pdf_list)):
        flux_ratio_dict[pdf_list[ii] + '_pdf'] = data_list[ii]

        error, xpeak = compute_onesig_pdf(data_list[ii], ratio_list[ii],
                                          usepeak=True, silent=True,
                                          verbose=True)
        
        fr_xpeak_dict[pdf_list[ii]+'_xpeak'] = xpeak
        fr_error_dict[pdf_list[ii]+'_lowhigh_error'] = error

    new_dict_list = [flux_ratio_dict, fr_xpeak_dict, fr_error_dict]
    # will have all the data and xpeaks for all histograms wanted
    histo_dict = OrderedDict()
    for bb in range(len(new_dict_list)):
        # dictionary = np.load(path of the dictionary)
        dic0 = new_dict_list[bb]
        # Added small dictionary to big dictionary
        histo_dict.update(dic0)

    today = date.today()
    pdf_file = 'Flux_Ratio histogram_plots_' + "%02i%02i" \
               % (today.month, today.day)+'.pdf'
    histogram(path, histo_dict, input_file, pdf_file, valid_file,
              table_key='bin_ID', sharex=sharex)
