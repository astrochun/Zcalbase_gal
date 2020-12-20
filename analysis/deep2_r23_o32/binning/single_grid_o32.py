
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from ..logging import log_stdout


def single_grid_o32(pdf_pages, npz_outfile, R23, O32, galinbin, log=None):
    """
    Purpose:
      This file holds the function to bin data for one dimensional binning along
      log(O32) Not used in current analysis

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

    # One_dimensional binning for O32
    sort0 = np.argsort(O32)
    y_sort0 = O32[sort0]

    # 200 galaxies per bin
    n_bins = np.int(len(O32)/galinbin)
    log.info(f"n_bins: {n_bins}")

    # Initializing Arrays for npz_outfile later
    N_arr0 = np.zeros((1, n_bins))
    T_arr = np.zeros((1, n_bins), dtype=object)

    # Bin starts and stops initializing
    bin_start = np.zeros(n_bins)
    bin_end = np.zeros(n_bins)

    # Sets the bins
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*galinbin]
        bin_end[ii] = y_sort0[(ii+1)*galinbin-1]
        if ii == max(n_bins):
            bin_end[ii] = max(y_sort0)

    # Organizes data into bins
    for oo in range(n_bins):
        idx_arr = np.where((O32 >= bin_start[oo]) & (O32 <= bin_end[oo]))[0]
        N_arr0[0, oo] += len(idx_arr)
        T_arr[0, oo] = idx_arr

    # Plotting
    fig, ax = plt.subplots()

    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x = np.log10(R23[finite0])
    y = np.log10(O32[finite0])

    h_lines = np.log10(bin_start)
    h_lines_end = np.log10(bin_end)
    ax.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*',
               alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for val in h_lines:
        plt.axhline(y=val, linewidth=0.3, color='k')
    for val_end in h_lines_end:
        plt.axhline(y=val_end, linewidth=0.3, color='g')

    fig.savefig(pp, format='pdf')
    log.info(f"Writing: {pdf_pages}")
    pp.close()

    R23_grid = [np.min(R23)]

    log.info(f"Writing: {npz_outfile}")
    np.savez(npz_outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=bin_start,
             N_arr0=N_arr0)

    fig.clear()

    log.info("finished ...")
