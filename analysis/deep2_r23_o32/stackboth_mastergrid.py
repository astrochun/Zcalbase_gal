"""
THIS FUNCTION DOES THE STACKING FOR THE ALL BINNING METHODS EXCEPT VORNOI

In this file, I define the stacking code in a function that runs over the master grid
Creates a pdf and a fits file
fits file for this code is given the name of the PDF file

Emission lines in spectrum (not all being used currently in study) See MSC for subset
[3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46,
4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]
"""

import numpy as np
import numpy.ma as ma

from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

from os.path import exists, join
from astropy.convolution import Box1DKernel, convolve

from . import general
from .logging import log_stdout

from Metallicity_Stack_Commons import lambda0
from Metallicity_Stack_Commons.column_names import filename_dict


def movingaverage_box1d(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth


def master_stacking(fitspath, fitspath_ini, dataset, grid_data_file, name,
                    mask=True, log=None):
    """

    Purpose
    Function stacks all spectra in a given bin and produces tables of properties of that bin

    Parameters
    :param fitspath: str. save location of the current run
    :param fitspath_ini: str. save location of all of Zcalbase
    :param dataset: str. keyword used to define binning method
    :param grid_data_file: str. npz file that holds the information from the binning process
    :param name: str. name of the outputted pdf file with graphs
    :param mask: bool. optional input used to mask the night sky lines if inputted (default: None)
    :param log: LogClass. Default use log_stdout()
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    RestframeMaster = join(fitspath_ini, 'DEEP2_Commons/Images/Master_Grid.fits')
    log.info(f"Reading: {RestframeMaster}")
    image2D, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])

    pdf_pages = PdfPages(join(fitspath, name))  # open pdf document

    individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, \
        data3 = general.get_det3(fitspath, fitspath_ini, log=log)

    log.info(f"Reading: {grid_data_file}")
    grid_data = np.load(grid_data_file, allow_pickle=True)  # This is the npz file
    R23_minimum = grid_data['R23_minimum']
    O32_minimum = grid_data['O32_minimum']

    image2DM = np.nan_to_num(image2D[det3])

    if mask:
        mask_file = join(fitspath_ini,
                         'stacking_masks/MastermaskArray.fits')
        log.info(f"Reading: {mask_file}")
        maskM = fits.getdata(mask_file)
        image2DM = np.ma.masked_array(image2DM, maskM[det3])
        log.debug("Mask[det3]", maskM[det3])

    outfile = join(fitspath, filename_dict['comp_spec'])
    if not exists(outfile):
        stack_2d = np.zeros((len(R23_minimum) * len(O32_minimum), len(wave)),
                            dtype=np.float64)
    else:
        log.info(f"Reading: {outfile}")
        stack_2d = fits.getdata(outfile)

    avg_R23 = np.zeros(len(R23_minimum) * len(O32_minimum))   # Same as xBar
    avg_O32 = np.zeros(len(R23_minimum) * len(O32_minimum))   # Same as yBar
    R23_node = np.zeros(len(R23_minimum) * len(O32_minimum))  # Same as R23_minimum
    O32_node = np.zeros(len(R23_minimum) * len(O32_minimum))  # Same as O32_minimum
    R23_med = np.zeros(len(R23_minimum) * len(O32_minimum))   # median R23 value
    O32_med = np.zeros(len(R23_minimum) * len(O32_minimum))   # median O32 value
    N_gal = np.zeros(len(R23_minimum) * len(O32_minimum))     # Same as Number_inbin

    n_N = R23_minimum.shape[0]
    if dataset in ['n_Bins', 'Double_Bin']:
        n_M = R23_minimum.shape[1]
    else:
        n_M = O32_minimum.shape[0]
    count = 0

    for rr in range(n_N):
        for oo in range(n_M):
            log.info(f"{rr}, {oo}")
            index = grid_data['locator'][rr, oo]

            # calculating the average and minimum values of R23 and O32
            if len(index) > 10:
                R23_node[count] = R23_minimum[rr, oo]
                O32_node[count] = O32_minimum[rr, oo]
                avg_R23[count] = np.average(R23[index])  # np.log10(R23)
                avg_O32[count] = np.average(O32[index])  # (np.log10(O32
                R23_med[count] = np.median(R23[index])   # np.log10(R23)
                O32_med[count] = np.median(O32[index])
                N_gal[count] = len(index)
                subgrid = image2DM[index]

                log.info(f"R23: {R23_node[count]} O32: {O32_node[count]}," +
                         f"avg_R23: {avg_R23[count]} avg_O32: {avg_O32[count]}")

                if exists(outfile):
                    Spect1D = stack_2d[count]
                else:
                    if mask:
                        Spect1D = np.ma.mean(subgrid, axis=0)
                    else:
                        Spect1D = np.nanmean(subgrid, axis=0)

                    stack_2d[count] = Spect1D

                # Compute number of spectra at a given wavelength
                a = ma.count(subgrid, axis=0)

                fig, ax = plt.subplots()

                # GridSpec
                gs = GridSpec(6, 1)
                ax1 = plt.subplot(gs[0, 0])
                ax2 = plt.subplot(gs[1:, 0])

                # ax1 Plot
                ax1.plot(wave, a, linewidth=0.5)

                # ax2 plot
                ax2.plot(wave, Spect1D, linewidth=0.5)
                ax2.set_xlabel('Wavelength')
                ax2.set_ylabel('Spectra 1D')
                ax2.minorticks_on()
                ax2.set_ylim([0, np.nanmax(Spect1D)])

                ax2.legend(loc='upper left', numpoints=3)

                str0 = f"R23={R23_minimum[rr, oo]:.1f} " + \
                       f"O32={O32_minimum[rr, oo]:.1f} " + \
                       f"N={len(index):d}"
                ax2.annotate(str0, (0.05, 0.95), xycoords='axes fraction',
                             ha='left', va='top', weight='bold')

                for x in lambda0:
                    ax2.axvline(x=x, linewidth=0.3, color='k')

                fig.set_size_inches(8, 11)

                x1, x2 = 4200, 4500
                ind = np.where((wave >= x1) & (wave <= x2))[0]
                sig, med = np.nanstd(Spect1D[ind]), np.nanmedian(Spect1D[ind])
                y1, y2 = med - 2.5 * sig, np.nanmax(Spect1D[ind])
                axins = zoomed_inset_axes(ax2, 2, loc=6,
                                          bbox_to_anchor=[375, 425])
                axins.plot(wave, Spect1D, linewidth=0.5)
                axins.set_xlim([x1, x2])
                axins.set_ylim([y1, y2])
                axins.minorticks_on()

                for x in lambda0:
                    axins.axvline(x=x, linewidth=0.3, color='r')

                mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="k",
                           ls='dashed', lw=0.5)

                fig.tight_layout()
                plt.draw()
                pdf_pages.savefig(fig)
                count += 1

    log.info(f"Writing: {join(fitspath, name)}")
    pdf_pages.close()

    # Writing fits file
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] = 1.00000
    header['CDELT2'] = 1.00000
    header['CRPIX2'] = 1.00000

    log.info(f"Writing: {outfile}")
    fits.writeto(outfile, stack_2d[0:count], header, overwrite=True)

    # Writing Ascii Tables and Fits Tables
    out_ascii = join(fitspath, filename_dict['bin_info'])  # Used to be 'binning_averages.tbl'

    ID = np.arange(0, len(R23_node), 1, dtype=int)
    n = ('bin_ID', 'logR23_min', 'logO32_min', 'logR23_avg', 'logO32_avg',
         'logR23_med', 'logO32_med', 'N_stack')
    # for n_split xnode and ynode are the lowest values of the bin while xBar and yBar are the averages

    tab0 = Table([ID, R23_node, O32_node, avg_R23, avg_O32, R23_med, O32_med, N_gal],
                 names=n)
    log.info(f"Writing: {out_ascii}")
    asc.write(tab0[0:count], out_ascii, format='fixed_width_two_line')

    fig.clear()

    log.debug("finished ...")
