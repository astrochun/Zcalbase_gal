from os.path import join, exists

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table

from ..log_commons import log_stdout


def two_times_binned(fitspath, pdf_file, npz_outfile, R23, O32, SNR3, data3,
                     galinbin, log=None):
    """
    One_dimensional binning for R23 followed by each bin being split in
    half one with high O32 and one with low O32

    Followed by plotting the spectra in their R23-O32 bins for visual

    :param fitspath: str. Path where files are retrieved and saved to
    :param pdf_file: str. Name of outputted pdf file
                     fitspath + name_dict['gridpdf_suffix']
    :param npz_outfile: str. Name of the npz file produced by the function
                        fitspath + name_dict['gridnpz_suffix']
    :param R23: np.array. Array of R23 from the get_det3 function
    :param O32: np.array. Array of O32 from the get_det3 function
    :param SNR3: np.array. Array of SNR3 from the get_det3 function
    :param data3: np.array. From get_det3 - indicates if spectra can be used
    :param galinbin: str. Array that specifies how many spectra in each bin
    :param log: LogClass or logging object

    No returns
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")
    pp = PdfPages(pdf_file)

    # One_dimensional binning for R23

    sort0 = np.argsort(R23)
    R23_sort0 = R23[sort0]

    n_bins = len(galinbin)
    
    log.info(f"n_bins: {n_bins}")
    n_bins_range = np.arange(0, 2*n_bins, 1)

    # Initializing Arrays for Grid stacking
    N_arr0 = np.zeros((n_bins, 2))
    T_arr = np.zeros((n_bins, 2), dtype=object)
    O32_grid = np.zeros((n_bins, 2))
    R23_grid = np.zeros((n_bins, 2))
    xBar = np.zeros(2 * n_bins)
    yBar = np.zeros(2 * n_bins)
    area = np.zeros(2 * n_bins)
    N_bin = np.zeros(len(data3), dtype=int)

    # Bin starts and stops initializing
    bin_start_1 = np.zeros(n_bins, dtype=np.int)
    bin_end_1 = np.zeros(n_bins, dtype=np.int)

    # Sets the bins increase the number of galaxies as R23 increases
    log.info(f"galinbin: {galinbin}")
    for ii in range(n_bins):
        if ii == 0:
            bin_start_1[ii] = 0
        else:
            bin_start_1[ii] = bin_end_1[ii-1]+1
        if ii == n_bins-1:
            bin_end_1[ii] = len(R23_sort0)-1
        else:
            log.info(galinbin[ii] + bin_start_1[ii])
            bin_end_1[ii] = galinbin[ii]+bin_start_1[ii]-1   
        log.info(f"Bin Start: {bin_start_1[ii]}  Bin end: {bin_end_1[ii]}")

        idx1 = np.where((R23 >= R23_sort0[bin_start_1[ii]]) &
                        (R23 <= R23_sort0[bin_end_1[ii]]))[0]
        med_val = np.median(O32[idx1])

        idx2 = np.where((R23 >= R23_sort0[bin_start_1[ii]]) &
                        (R23 <= R23_sort0[bin_end_1[ii]]) &
                        (O32 <= med_val))[0]
        idx3 = np.where((R23 >= R23_sort0[bin_start_1[ii]]) &
                        (R23 <= R23_sort0[bin_end_1[ii]]) &
                        (O32 > med_val))[0]

        O32_grid[ii, 0] = np.median(O32[idx2])
        O32_grid[ii, 1] = np.median(O32[idx3])
        R23_grid[ii, 0] = R23_sort0[bin_start_1[ii]]
        R23_grid[ii, 1] = R23_sort0[bin_start_1[ii]]

        N_arr0[ii, 0] += len(idx2)
        N_arr0[ii, 1] += len(idx3)
        
        T_arr[ii, 0] = idx2
        T_arr[ii, 1] = idx3
        
        N_bin[idx2] = ii*2
        N_bin[idx3] = ii*2 + 1

        xBar[ii*2] = np.log10(np.average(R23[idx2]))
        xBar[ii*2+1] = np.log10(np.average(R23[idx3]))

        yBar[ii*2] = np.log10(np.average(O32[idx2]))
        yBar[ii*2+1] = np.log10(np.average(O32[idx3]))
        
        area[ii*2] = len(idx2)
        area[ii*2+1] = len(idx3)

    O32_values = O32_grid.reshape(2 * n_bins)
    R23_values = R23_grid.reshape(2 * n_bins)
    log.info(f"O32_values: {O32_values}")
    log.info(f"R23_values: {R23_values}")
        
    # Plotting
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    h_lines = np.log10(R23_values)
    ax.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for pp in range(2*len(bin_start_1)):
        plt.axvline(x=h_lines[pp], linewidth=0.3, color='k')
    fig.savefig(pp, format='pdf')

    log.info(f"Writing: {pdf_file}")
    pp.close()
    fig.clear()

    # Save NPZ
    log.info(f"Writing: {npz_outfile}")
    np.savez(npz_outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid,
             N_arr0=N_arr0)

    n1 = ('ID', 'R23_value', 'O32_median', 'xBar', 'yBar', 'area')
    tab1 = Table([n_bins_range, R23_values, O32_values, xBar, yBar, area],
                 names=n1)
    out_table_file = join(fitspath, "Double_Bin_binning_averages.tbl")
    if not exists(out_table_file):
        log.info(f"Writing: {out_table_file}")
    else:
        log.info(f"Overwriting : {out_table_file}")
    asc.write(tab1, out_table_file, format='fixed_width_two_line', overwrite=True)

    n2 = ('R23', 'O32', 'SN_5007', 'N_bin')
    tab2 = Table([R23, O32, SNR3, N_bin], names=n2)
    out_table_file3 = join(fitspath, "Double_Bin_2d_binning_datadet3.tbl")
    if not exists(out_table_file3):
        log.info(f"Writing: {out_table_file3}")
    else:
        log.info(f"Overwriting : {out_table_file3}")
    asc.write(tab2, out_table_file3, format='fixed_width_two_line', overwrite=True)

    log.debug("finished.")
