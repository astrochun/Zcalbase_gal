import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from os.path import join, exists

from ..plotting import manual_bin_grid

from Metallicity_Stack_Commons.column_names import filename_dict

from ..log_commons import log_stdout


def n_times_binned(fitspath, pdf_pages, npz_outfile, n_split, individual_ID,
                   R23, O32, SNR3, data3, galinbin, log=None):
    """
    This file holds the function to bin data adaptively based on the entered
    number of galaxies for each bin.
    First, it bins in the R23 direction and then in the O32 direction

    This is used in current analysis.

    :param fitspath: str. Absolute path for working directory
    :param pdf_pages: str. name of outputted pdf file
                      fitspath + name_dict['gridpdf_suffix']
    :param npz_outfile: str. name of the npz file produced by the function
                        fitspath + name_dict['gridnpz_suffix']
    :param n_split: int. the number of times each R23 bin will be split
    :param galinbin: np.array. Number of spectra in each bin
    :param individual_ID: np.array. Individual IDs for each spectra
    :param R23: np.array. log(R23) measurements of each spectra
    :param O32: np.array. log(O32) measurements of each spectra
    :param SNR3: np.array. Signal to noise of the OIII emission lines
                            of each spectra
    :param data3: np.array. From get_det3 - indicates if spectra can be used
    :param log: LogClass or logging object

    No returns
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    pp = PdfPages(pdf_pages)

    # One_dimensional binning for R23
    sortR23 = np.argsort(R23)
    R23_sort0 = R23[sortR23]

    n_bins = len(galinbin)
    n_bins_range = np.arange(0, n_split*n_bins, 1)

    # Initializing Arrays for Grid stacking
    inbin_shape = (n_bins, n_split)
    n_bins_total = n_split * n_bins
    number_inbin = np.zeros(inbin_shape)                 # Used to be N_arr0
    locator = np.zeros((n_bins, n_split), dtype=object)  # Used to be T_arr
    O32_minimum = np.zeros(inbin_shape)                  # Used to be O32_grid
    R23_minimum = np.zeros(inbin_shape)                  # Used to be R23_grid
    O32_median = np.zeros(inbin_shape)
    R23_median = np.zeros(inbin_shape)
    O32_max = np.zeros(inbin_shape)
    R23_max = np.zeros(inbin_shape)
    xBar = np.zeros(n_bins_total)
    yBar = np.zeros(n_bins_total)
    area = np.zeros(n_bins_total)
    bin_number = np.zeros(len(data3), dtype=int)  # Used to be N_bin

    R23_minall = np.zeros(len(R23))
    O32_minall = np.zeros(len(R23))
    R23_avgall = np.zeros(len(R23))
    O32_avgall = np.zeros(len(R23))
    R23_medall = np.zeros(len(R23))
    O32_medall = np.zeros(len(R23))
    R23_maxall = np.zeros(len(R23))
    O32_maxall = np.zeros(len(R23))

    # Bin starts and stops initializing
    bin_start = np.zeros(n_bins, dtype=np.int)
    bin_end = np.zeros(n_bins, dtype=np.int)

    for ii in range(n_bins):
        log.info(f"ii: {ii}")
        if ii == 0:
            bin_start[ii] = 0
        else:
            bin_start[ii] = bin_end[ii-1] + 1
        if ii == n_bins - 1:
            bin_end[ii] = len(R23_sort0) - 1
        else:
            bin_end[ii] = galinbin[ii]+bin_start[ii] - 1
        log.info(f"Bin Start: {bin_start[ii]}  Bin end: {bin_end[ii]}")

        R23_idx = np.where((R23 >= R23_sort0[bin_start[ii]]) &
                           (R23 <= R23_sort0[bin_end[ii]]))[0]

        # O32 relative to the R23_index/ O32 sorted into the R23 bins
        O32_inbins = O32[R23_idx]
        # Sort the O32 in their bins so that we have their locations
        O32_index = np.argsort(O32_inbins)
        # Take O32 in each bin & organizes them based on index in previous line
        sortO32 = O32_inbins[O32_index]

        # The following lines could go into a definition called
        # "evenSplit( ... )" therefore allowing you to create a different
        # splitting method definition called "optimalSplit( ... )"
        # or something to that effect if so desired.
        n_subbins = np.int(np.floor(float(len(sortO32))/n_split))
        subbin_arr = np.ones(n_split, dtype=int) * n_subbins
        n_remainder = len(sortO32) - n_subbins*n_split
        backwardsIdx = -1

        # Sorting all the remainders into the correct bins:
        # Use backwardsIdx or forwardsIdx to start adding spectra
        # from the front or the back of the list of bins
        for rr in range(n_remainder):
            subbin_arr[backwardsIdx] = n_subbins + 1
            backwardsIdx -= 1
        subbin_arr[backwardsIdx] = n_subbins
        
        startIdx = 0
        endIdx = subbin_arr[0]
        for jj in range(n_split):
            # Let's grab all O32 values
            O32_values_perbin = sortO32[startIdx:endIdx]
            O32_inbins_idx = O32_index[startIdx:endIdx]

            # This index gives the positions of all the R23 values relative
            # to the O32 values over the entire R23 bin
            N_bin_idx = R23_idx[O32_inbins_idx]        

            # Gives the bin number for each spectra
            bin_number[N_bin_idx] = (ii * n_split) + jj

            # Now let's start sorting our data into variables to use later

            # First two map minimum R23 and O32 measure
            # to all individual spectra in bin
            # Second two values give the lowest O32 and R23 value set for each bin
            R23_minall[N_bin_idx] = R23_sort0[bin_start[ii]]
            O32_minall[N_bin_idx] = O32_values_perbin[0]
            R23_minimum[ii, jj] = R23_sort0[bin_start[ii]]
            O32_minimum[ii, jj] = O32_values_perbin[0]

            # First two map median R23 and O32 measure
            # to all individual spectra in bin
            # Second two give the median R23 and O32 value for each bin
            R23_medall[N_bin_idx] = np.median(R23[N_bin_idx])
            O32_medall[N_bin_idx] = np.median(O32_values_perbin)
            R23_median[ii, jj] = np.median(R23[N_bin_idx])
            O32_median[ii, jj] = np.median(O32_values_perbin)

            # First two map maximum R23 and O32 measure
            # to all individual spectra in bin
            # Second two give the maximum R23 and O32 value for each bin
            R23_maxall[N_bin_idx] = np.max(R23[N_bin_idx])
            O32_maxall[N_bin_idx] = np.max(O32[N_bin_idx])
            R23_max[ii, jj] = np.max(R23[N_bin_idx])
            O32_max[ii, jj] = np.max(O32[N_bin_idx])

            # Gives the number of galaxies in each bin
            number_inbin[ii, jj] += len(O32_values_perbin)

            # Gives the index (location numbers) for the spectra in each bin
            # and is used later loop over and get the galaxies
            locator[ii, jj] = N_bin_idx

            # Maps average R23 measure to all individual spectra in bin
            R23_avgall[N_bin_idx] = np.average(R23[N_bin_idx])
            # Gives the average R23 value for each bin
            xBar[(ii*n_split)+jj] = np.average(R23[N_bin_idx])

            # Maps average O32 measure to all individual spectra in bin
            O32_avgall[N_bin_idx] = np.average(O32_values_perbin)
            # Gives the average O32 value for each bin
            yBar[(ii*n_split)+jj] = np.average(O32_values_perbin)

            # Gives the number of galaxies in each bin
            area[(ii*n_split)+jj] = len(O32_values_perbin)

            # Now we can shift our window over for the next split
            # like the counting index the Stacking code
            startIdx = endIdx
            endIdx = startIdx+subbin_arr[jj] + 1

    O32_lowlimit = O32_minimum.reshape(n_bins_total)  # Used to be O32_grids
    R23_lowlimit = R23_minimum.reshape(n_bins_total)  # Used to be R23_grids
    O32_medians = O32_median.reshape(n_bins_total)
    R23_medians = R23_median.reshape(n_bins_total)
    R23_maxval = R23_max.reshape(n_bins_total)
    O32_maxval = O32_max.reshape(n_bins_total)

    # Plotting
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    vlines = np.log10(R23_lowlimit)
    hlines = np.log10(O32_lowlimit)
    ax.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')

    for jj in range(len(O32_lowlimit)):
        x_min = vlines[jj]
        if jj <= (len(O32_lowlimit) - n_split - 1):
            x_max = vlines[jj + n_split]
        else:
            x_max = np.log10(max(R23))
        plt.axvline(x=vlines[jj], linewidth=0.3, color='k')
        
        x_value = [x_min, x_max]
        y_value = [hlines[jj], hlines[jj]]
        y_average = [yBar[jj], yBar[jj]]
        plt.plot(x_value, y_value, linewidth=0.3, color='b')
        plt.plot(x_value, y_average, linewidth=0.3, color='g')

    fig.savefig(pp, format='pdf')
    log.info(f"Writing: {pdf_pages}")
    pp.close()

    # New plotting function
    manual_bin_grid.graph_bins(fitspath, n_split)

    # Writing as a dictionary
    avg_dict = {'locator': locator,
                'R23_minimum': R23_minimum,
                'O32_minimum': O32_minimum,
                'Number_inbin': number_inbin}
    log.info(f"Writing: {npz_outfile}")
    np.savez(npz_outfile, **avg_dict)

    n1 = ('bin_ID', 'N_stack', 'logR23_min', 'logO32_min', 'logR23_avg',
          'logO32_avg', 'logR23_median', 'logO32_median', 'logR23_max',
          'logO32_max')
    tab1 = Table([n_bins_range, area, np.log10(R23_lowlimit),
                  np.log10(O32_lowlimit), np.log10(xBar), np.log10(yBar),
                  np.log10(R23_medians), np.log10(O32_medians),
                  np.log10(R23_maxval), np.log10(O32_maxval)], names=n1)

    # Used to be called +dataset+'_binning_averages.tbl
    bin_outfile = join(fitspath, filename_dict['bin_info'])
    if not exists(bin_outfile):
        log.info(f"Writing: {bin_outfile}")
    else:
        log.info(f"Overwriting : {bin_outfile}")
    asc.write(tab1, bin_outfile, format='fixed_width_two_line', overwrite=True)

    n2 = ('logR23', 'logO32', 'OIII_5007_S/N', 'bin_ID', 'ID', 'logR23_min',
          'logO32_min', 'logR23_avg', 'logO32_avg',
          'logR23_median', 'logO32_median', 'logR23_max', 'logO32_max')
    tab2 = Table([R23, O32, SNR3, bin_number, individual_ID,
                 R23_minall, O32_minall, R23_avgall, O32_avgall,
                 R23_medall, O32_medall, R23_maxall, O32_maxall], names=n2)

    # Used to be + dataset+'_2d_binning_datadet3.tbl
    indv_outfile = join(fitspath, filename_dict['indv_bin_info'])
    if not exists(indv_outfile):
        log.info(f"Writing: {indv_outfile}")
    else:
        log.info(f"Overwriting : {indv_outfile}")
    asc.write(tab2, indv_outfile, format='fixed_width_two_line', overwrite=True)

    '''n3 = ('ID' , 'R23_grid', 'O32_grid')
    tab1 = Table([n_bins_range, R23_grid, O32_grid], names = n3)
    log.info(f"Writing: {fitspath+'/Double_Bin_grid_values.tbl'}")
    log.info(f"Overwriting : {fitspath+'/Double_Bin_grid_values.tbl'}")
    asc.write(tab1, fitspath+'/Double_Bin_grid_values.tbl', format='fixed_width_two_line', overwrite=True)'''
    # Create another ascii table with the R23_grid and O32_grid values for plots

    log.debug("finished.")
