# THIS CODE IS NOT USED IN THE GENERAL FUNCTION AND IS NOT CURRENTLY BEING RUN

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust


def redshift_graphs(fitspath, RestframeMaster, binning_grid_npz, DEEP2_all_files):
    """
    Purpose
    This code makes the multi graph and then graphs the redshifts in each bin
    Used as an error checker to make sure results were reasonable

    Parameters
    ----------
    fitspath ->  path where files are called from and saved to
    RestframeMaster -> rest frame for all spectra
                    ie. RestframeMaster = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
    binning_grid_npz -> npz file produced by binning code
    DEEP2_all_files -> raw data files
    """
    pp = PdfPages(fitspath + 'Redshift_Binning.pdf')
    pp2 = PdfPages(fitspath + 'OH_masking_check.pdf')

    maskM, maskHdr = fits.getdata(RestframeMaster, header=True)  # update header

    grid_data = np.load(binning_grid_npz)

    tab, header = fits.getdata(DEEP2_all_files, header=True)

    R23_grid = grid_data['R23_minimum']
    O32_grid = grid_data['O32_minimum']

    nrows = 4
    ncols = 4

    xlim = [np.min(tab.ZSPEC), np.max(tab.ZSPEC)]

    count = 0

    wave0 = maskHdr['CRVAL1'] + maskHdr['CDELT1'] * np.arange(maskHdr['NAXIS1'])

    for rr in range(len(R23_grid)):
        for oo in range(len(O32_grid)):
    
            index = grid_data['T_arr'][rr, oo]
            if len(index) > 10:
                if count % (nrows*ncols) == 0:
                    fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols)
                    # endif
                
                row = count / nrows % ncols
                col = count % ncols
        
                t_ax = ax_arr[row, col]
                t_zspec = tab.ZSPEC[index]

                # Graphing
                hist0 = t_ax.hist(t_zspec, bins=30)

                txt0 = r'R$_{23}$=%.3f  O$_{32}$=%.3f  N=%i' % (R23_grid[rr], O32_grid[oo], len(t_zspec)) + '\n'
                txt0 += 'Avg: %.3f  Med: %.3f\n' % (np.average(t_zspec), np.median(t_zspec))
                txt0 += 'Min: %.3f  Max: %.3f\n' % (np.min(t_zspec), np.max(t_zspec))
                txt0 += r'$\sigma$: %.3f' % np.std(t_zspec)
                t_ax.annotate(txt0, [0.95, 0.95], xycoords='axes fraction', va='top', ha='right')
                # t_ax.set_title(str('R23='+str(R23_grid[rr]))+ '&'+ 'O32='+str(O32_grid[oo]))
                t_ax.set_xlim(xlim)
                t_ax.set_ylim([0, np.max(hist0[0]) * 1.25])
                # t_ax.set_ylim([0,70])
            
                if row != nrows-1:
                    t_ax.set_xticklabels([])
                else:
                    t_ax.set_xlabel('Redshift')
                    # ax_arr.ylabel('y-axis')
                
                if count % (nrows * ncols) == nrows * ncols - 1:
                    print(count)
                    subplots_adjust(left=0.2, right=0.98, bottom=0.06, top=0.95,
                                    hspace=0.0)
                    fig.set_size_inches(11, 8)
                    plt.draw()
                    fig.savefig(pp, format='pdf')
                    # pdf_file.savefig(fig)

                if count % nrows == 0:
                    fig2, ax_arr2 = plt.subplots(nrows=nrows, ncols=1)

                OHmask_arr = np.zeros(len(wave0))
            
                for tt in xrange(len(t_zspec)):
                    msk = np.where(maskM[index[tt]] == 1)[0]
                    OHmask_arr[msk] += 1
                
                row = count % nrows
                ax_arr2[row].plot(wave0, len(index)-OHmask_arr)
                if count % nrows == nrows-1:
                    print('blah2: ', count)
                    subplots_adjust(left=0.03, right=0.98, bottom=0.06, top=0.95, hspace=0.0)
                    fig2.set_size_inches(11, 8)
                    plt.draw()
                    fig2.savefig(pp2, format='pdf')

                count += 1

    for cc in range(count % (nrows*ncols), nrows*ncols):
        row0 = cc / nrows % ncols
        col0 = cc % ncols
        ax_arr[row0, col0].axis('off')
        # ax_arr[row0,col0].set_xticklabels()

    if count % (nrows * ncols) != nrows * ncols - 1:
        for cc in range(ncols):
            ax_arr[row, cc].set_xlabel('Redshift')
            fig.set_size_inches(11, 8)
            subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
            pp.savefig(fig)

        fig2.set_size_inches(11, 8)
        subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
        pp2.savefig(fig2)

    pp.close()
    pp2.close()
