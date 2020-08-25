
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def single_grid_r23(fitspath, pdf_pages, outfile, R23, O32, galinbin):
    """
    This file holds the function to bin data for one dimensional binning along R23
    Not used in current analysis

    Inputs:
    fitspath  -> path where files are called from and saved to
    pdf_pages -> name of outputted pdf file
    outfile   -> name of the npz file produced by the function
    galinbin  -> array of numbers that specifies how many spectra go in each bin
    Other variables -> emission file values of spectra that come from the get_det3 function
    """
    pdf_pages = PdfPages(pdf_pages)

    # One_dimensional binning for R23

    sort0 = np.argsort(R23)
    R23_sort0 = R23[sort0]

    # 200 galaxies per bin
    n_bins = len(galinbin)
    print n_bins

    # Initializing Arrays for outfile later
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
            print galinbin[ii]+bin_start[ii]
            bin_end[ii] = galinbin[ii]+bin_start[ii]-1   
        print 'Bin Start:', bin_start[ii], 'Bin end:', bin_end[ii]

        idx_arr = np.where((R23 >= R23_sort0[bin_start[ii]]) & (R23 <= R23_sort0[bin_end[ii]]))[0]
        R23_grid[ii, 0] = R23_sort0[bin_start[ii]]
        N_arr0[ii, 0] += len(idx_arr)
        T_arr[ii, 0] = idx_arr
        print('Start:', R23_grid[ii])   # This could be incorrect so I need to check this again###

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
        plt.axvline(x=hlines[pp], linewidth=0.3, color='k')
    for ll in range(len(bin_end)):
        plt.axvline(x=hlines_end[ll], linewidth=0.3, color='g')
    fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()

    O32_grid = [np.min(O32)]
    print('R23_grid:', bin_start)
    print('O32_grid:', O32_grid)
    print('T_arr:', T_arr)
    print('N_arr0:', N_arr0)

    np.savez(outfile, T_arr=T_arr, R23_grid=bin_start, O32_grid=O32_grid, N_arr0=N_arr0)

    fig.clear()
