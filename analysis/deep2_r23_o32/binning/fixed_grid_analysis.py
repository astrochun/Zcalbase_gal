
# Currently binning with bin=0.25

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from astropy.table import Table
from os.path import join


def making_grid(fitspath, pdf_pages, outfile, R23, O32, det3, R23_bin, O32_bin):
    """
    This file holds the function to bin data in a bin of fixed size entered as an input
    Not used in current analysis

    Inputs:
    fitspath  -> path where files are called from and saved to
    pdf_pages -> name of outputted pdf file
    outfile   -> name of the npz file produced by the function
    R23_bin   -> set size of the bin in the R23 direction
    O32_bin   -> set size of the bin in the O32 direction
    Other variables -> emission file values of spectra that come from the get_det3 function
    """

    fig1, ax1 = plt.subplots() 

    xlim = [0.4, 50]
    ylim = [0.1, 20]

    R23_grid = np.arange(np.log10(xlim[0]), np.log10(xlim[1])+R23_bin, R23_bin)
    O32_grid = np.arange(np.log10(ylim[0]), np.log10(ylim[1])+R23_bin, O32_bin)
    N_arr0 = np.zeros((len(R23_grid), len(O32_grid)), dtype=np.float)

    N = len(R23_grid)
    M = len(O32_grid)
    T_arr = np.zeros((N, M), dtype=object)

    # Plotting
    label0 = 'Field = data0[det3], N='+str(len(det3))
    x = np.log10(R23)
    y = np.log10(O32)
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x = x[finite0]
    y = y[finite0]
    ax1.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1, label=label0)
    x0 = x.tolist()
    y0 = y.tolist()

    for jj in range(len(R23_grid)):
        for kk in range(len(O32_grid)):
            array = np.where((x < R23_grid[jj]+R23_bin) & (x >= R23_grid[jj]) &
                             (y < O32_grid[kk]+O32_bin) & (y >= O32_grid[kk]))[0]
            print(array)
            N_arr0[jj, kk] += len(array)
            T_arr[jj, kk] = array

    ax1.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax1.set_xlabel(r'log($R_{23}$)')
    ax1.set_ylabel(r'log($O_{32}$)')
    ax1.set_xlim(np.log10(xlim))
    ax1.set_ylim(np.log10(ylim))
    ax1.minorticks_on()
    ax1.legend(loc='upper left', numpoints=3)     
    fig1.set_size_inches(9, 9)

    pdf_pages.savefig(fig1)

    fig2 = plt.figure()
    ax2 = plt.gca()  # ax = plt.subplots() #plt.gcf()
    cm = plt.cm.get_cmap('Blues')

    # Colorbar and hexbin plotting
    tabmastergrid = Table([x0, y0])
    tabmastergrid_name = join(fitspath, 'testmastergrid.tbl')
    asc.write(tabmastergrid, tabmastergrid_name, format='fixed_width_two_line')
    hb0 = ax2.hexbin(x0, y0, gridsize=(len(R23_grid), len(O32_grid)), cmap=cm)
    cbaxes = fig2.add_axes([0.135, 0.20, 0.75, 0.025])
    tick = [0., 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]
    cb = fig2.colorbar(hb0, cax=cbaxes, ticks=tick, orientation='horizontal')
    cb.set_label('density')

    ax2.scatter(x, y, 1.5, facecolor='r', edgecolor='face', marker='*', alpha=1, label=label0)
    ax2.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax2.set_xlabel(r'log($R_{23}$)')
    ax2.set_ylabel(r'log($O_{32}$)')
    ax2.set_xlim(np.log10(xlim))
    ax2.set_ylim(np.log10(ylim))
    ax2.minorticks_on()
    ax2.legend(loc='upper left', numpoints=3)

    fig2.set_size_inches(8, 8)
    pdf_pages.savefig(fig2)

    pdf_pages.close()

    fig1.clear()
    fig2.clear()

    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid, N_arr0=N_arr0)
