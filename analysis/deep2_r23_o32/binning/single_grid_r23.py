
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from ..log_commons import log_stdout


def single_grid_r23(pdf_pages, npz_outfile, R23, O32, galinbin, log=None):
    """
    Purpose:
      This file holds the function to bin data for one dimensional binning
      along R23.
      Note: Not used in current analysis

    :param pdf_pages: name of outputted pdf file
    :param npz_outfile: name of the npz file produced by the function
    :param R23: np.array. R23 measurements
    :param O32: np.array. O32 measurements
    :param galinbin: np.array. Number of spectra in each bin
    :param log: LogClass or logging object
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    pp = PdfPages(pdf_pages)

    # One_dimensional binning for R23

    sort0 = np.argsort(R23)
    R23_sort0 = R23[sort0]

    # 200 galaxies per bin
    n_bins = len(galinbin)
    log.info(f"n_bins: {n_bins}")

    # Initializing Arrays for npz_outfile later
    R23_grid = np.zeros((n_bins, 1))
    N_arr0 = np.zeros((n_bins, 1))
    T_arr = np.zeros((n_bins, 1), dtype=object)

    # Bin starts and stops initializing
    bin_start = np.zeros(n_bins, dtype=np.int)
    bin_end = np.zeros(n_bins, dtype=np.int)

    for ii in range(n_bins):
        if ii == 0:
            bin_start[ii] = 0
        else:
            bin_start[ii] = bin_end[ii-1]+1

        if ii == n_bins-1:
            bin_end[ii] = len(R23_sort0)-1
        else:
            bin_end[ii] = galinbin[ii]+bin_start[ii]-1

        log.info(f"Bin Start: {bin_start[ii]}  Bin end: {bin_end[ii]}")

        idx_arr = np.where((R23 >= R23_sort0[bin_start[ii]]) &
                           (R23 <= R23_sort0[bin_end[ii]]))[0]
        R23_grid[ii, 0] = R23_sort0[bin_start[ii]]
        N_arr0[ii, 0] += len(idx_arr)
        T_arr[ii, 0] = idx_arr
        log.info(f"Start: {R23_grid[ii]}")  # This could be incorrect so I need to check this again

    # Plotting
    fig, ax = plt.subplots()

    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x = np.log10(R23[finite0])
    y = np.log10(O32[finite0])

    h_lines = np.log10(bin_start)
    h_lines_end = np.log10(bin_end)
    ax.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')

    for val in h_lines:
        plt.axvline(x=val, linewidth=0.3, color='k')
    for val_end in h_lines_end:
        plt.axvline(x=val_end, linewidth=0.3, color='g')

    fig.savefig(pp, format='pdf')
    log.info(f"Writing: {pdf_pages}")
    pp.close()

    O32_grid = [np.min(O32)]
    log.info(f"R23_grid: {bin_start}")
    log.info(f"O32_grid: {O32_grid}")
    log.info(f"T_arr: {T_arr}")
    log.info(f"N_arr0: {N_arr0}")

    log.info(f"Writing: {npz_outfile}")
    np.savez(npz_outfile, T_arr=T_arr, R23_grid=bin_start, O32_grid=O32_grid,
             N_arr0=N_arr0)

    fig.clear()

    log.info("finished.")
