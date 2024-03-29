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

from . import name_dict, read_fitsfiles, get_det3
from .log_commons import log_stdout

from Metallicity_Stack_Commons import lambda0
from Metallicity_Stack_Commons.column_names import filename_dict


def master_stacking(fitspath, fitspath_ini, dataset, grid_data_file,
                    mask=True, log=None):
    """
    Function stacks all spectra in a given bin and produces tables of
    properties of that bin
    This function does the stacking for all binning methods except Vornoi

    :param fitspath: str. save location of the current run
    :param fitspath_ini: str. save location of all of Zcalbase
    :param dataset: str. keyword used to define binning method
    :param grid_data_file: str. npz file that holds the information
                           from the binning process
    :param mask: bool. optional input used to mask the night sky lines
                if inputted (default: None)
    :param log: LogClass. Default use log_stdout()

    PDF File: fitspath + stack_name
    TBL File: fitspath + filename_dict['bin_info']
    FITS File: fitspath + filename_dict['comp_spec']

    No returns
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    if mask:
        stack_name = name_dict['Stackname']
    else:
        stack_name = name_dict['Stackname_nomask']

    restframe_master_file = join(fitspath_ini,
                                 'DEEP2_Commons/Images/Master_Grid.fits')
    log.info(f"Reading: {restframe_master_file}")
    fits_dict = read_fitsfiles(restframe_master_file)
    image2D = fits_dict['fits_data']
    header = fits_dict['header']
    wave = fits_dict['wave']

    pdf_file = join(fitspath, stack_name)
    pp = PdfPages(pdf_file)  # open pdf document

    individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, \
        data3 = get_det3(fitspath, fitspath_ini, log=log)

    log.info(f"Reading: {grid_data_file}")
    # This is the npz file
    grid_data = np.load(grid_data_file, allow_pickle=True)
    R23_minimum = grid_data['R23_minimum']
    O32_minimum = grid_data['O32_minimum']
    len_R23 = len(R23_minimum)
    len_O32 = len(O32_minimum)

    image2DM = np.nan_to_num(image2D[det3])

    if mask:
        mask_file = join(fitspath_ini,
                         'DEEP2_Commons/Images/MastermaskArray.fits')
        log.info(f"Reading: {mask_file}")
        maskM = fits.getdata(mask_file)
        image2DM = np.ma.masked_array(image2DM, maskM[det3])
        log.debug("Mask[det3]", maskM[det3])

    stack_2d_file = join(fitspath, filename_dict['comp_spec'])
    if not exists(stack_2d_file):
        stack_2d_tab = np.zeros((len_R23 * len_O32, len(wave)),
                                dtype=np.float64)
    else:
        log.info(f"Reading: {stack_2d_file}")
        stack_2d_tab = fits.getdata(stack_2d_file)

    # Same as xBar
    avg_R23 = np.zeros(len(R23_minimum) * len(O32_minimum))
    # Same as yBar
    avg_O32 = np.zeros(len(R23_minimum) * len(O32_minimum))
    # Same as R23_minimum
    R23_node = np.zeros(len(R23_minimum) * len(O32_minimum))
    # Same as O32_minimum
    O32_node = np.zeros(len(R23_minimum) * len(O32_minimum))
    # median R23 value
    R23_med = np.zeros(len(R23_minimum) * len(O32_minimum))
    # median O32 value
    O32_med = np.zeros(len(R23_minimum) * len(O32_minimum))
    # Same as Number_inbin
    N_gal = np.zeros(len(R23_minimum) * len(O32_minimum))

    n_N = R23_minimum.shape[0]
    if dataset in ['n_Bins', 'Double_Bin']:
        n_M = R23_minimum.shape[1]
    else:
        n_M = O32_minimum.shape[0]
    count = 0

    for rr in range(n_N):
        for oo in range(n_M):
            log.info(f"rr: {rr}  oo: {oo}")
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

                log.info(f"R23: {R23_node[count]} O32: {O32_node[count]} " +
                         f"avg_R23: {avg_R23[count]} "
                         f"avg_O32: {avg_O32[count]}")

                if exists(stack_2d_file):
                    Spect1D = stack_2d_tab[count]
                else:
                    if mask:
                        Spect1D = np.ma.mean(subgrid, axis=0)
                    else:
                        Spect1D = np.nanmean(subgrid, axis=0)

                    stack_2d_tab[count] = Spect1D

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
                pp.savefig(fig)
                count += 1

    log.info(f"Writing: {pdf_file}")
    pp.close()

    # Writing fits file
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] = 1.00000
    header['CDELT2'] = 1.00000
    header['CRPIX2'] = 1.00000

    if not exists(stack_2d_file):
        log.info(f"Writing: {stack_2d_file}")
        fits.writeto(stack_2d_file, stack_2d_tab[0:count], header,
                     overwrite=True)

    # Writing Ascii Tables and Fits Tables
    # Used to be 'binning_averages.tbl'
    bin_info_file = join(fitspath, filename_dict['bin_info'])

    ID = np.arange(0, len(R23_node), 1, dtype=int)
    n = ('bin_ID', 'logR23_min', 'logO32_min', 'logR23_avg', 'logO32_avg',
         'logR23_med', 'logO32_med', 'N_stack')
    # for n_split xnode and ynode are the lowest values
    # of the bin while xBar and yBar are the averages

    bin_info_tab = Table([ID, R23_node, O32_node, avg_R23,
                          avg_O32, R23_med, O32_med, N_gal], names=n)
    if not exists(bin_info_file):
        log.info(f"Writing: {bin_info_file}")
    else:
        log.info(f"Overwriting: {bin_info_file}")
    asc.write(bin_info_tab[0:count], bin_info_file,
              format='fixed_width_two_line')

    fig.clear()

    log.debug("finished.")
