
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def single_grid_o32(fitspath, pdf_pages, outfile, R23, O32, galinbin):
    """
    This file holds the function to bin data for one dimensional binning along O32
    Not used in current analysis

    Inputs:
    fitspath  -> path where files are called from and saved to
    pdf_pages -> name of outputted pdf file
    outfile   -> name of the npz file produced by the function
    galinbin  -> array of numbers that specifies how many spectra go in each bin
    Other variables -> emission file values of spectra that come from the get_det3 function
    """
    pdf_pages = PdfPages(pdf_pages)

    # One_dimensional binning for O32

    sort0 = np.argsort(O32)
    y_sort0 = O32[sort0]

    # 200 galaxies per bin
    n_bins = np.int(len(O32)/galinbin)
    print n_bins

    # Initializing Arrays for outfile later
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
    x = R23
    y = O32
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x1 = x[finite0]
    y1 = y[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    hlines = np.log10(bin_start)
    hlines_end = np.log10(bin_end)
    ax.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for pp in range(len(bin_start)):
        plt.axhline(y=hlines[pp], linewidth=0.3, color='k')
    for ll in range(len(bin_end)):
        plt.axhline(y=hlines_end[ll], linewidth=0.3, color='g')
    fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()

    R23_grid = [np.min(R23)]
    
    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=bin_start, N_arr0=N_arr0)

    fig.clear()
