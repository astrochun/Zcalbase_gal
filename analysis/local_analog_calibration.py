from astropy.io import ascii as asc
import numpy as np
import matplotlib.pyplot as plt

from .deep2_r23_o32.log_commons import log_stdout, LogClass
from .deep2_r23_o32 import bian_coeff, ctype, bian18_R23_OH, bian18_O32_OH, \
    bian18_OH_O32
from .green_peas_calibration import get_zcalbase_sample
from .deep2_r23_o32.plotting import plot_difference_curvefit


def main(lR23, lO32, OH, out_pdf, R23_pdf_file, n_bins=4, ID=[],
         lR23_err=[], lO32_err=[], OH_err=[], R23_xra=[], O32_xra=[],
         yra=[], ctype=[], label=[''], marker=[], log=None):
    """
    Plots R23 and O32 line ratios and compare it to Bian et al. (2018)
    calibration that is based on local analogs to z~2 galaxies
    Main function to plot dataset against Bian+ (2018) calibration

    :param lR23: array. log(R_23)
    :param lO32: array. log([OIII]/[OII])
    :param OH: array. 12+log(O/H)
    :param out_pdf: str. full path for output PDF
    :param ID: array. list of list of ID's
    :param lR23_err: array. error on x-axis (logR23)
    :param lO32_err: array. error on x-axis (logO32)
    :param OH_err: array. metallicity error on y-axis
    :param R23_xra: array. x-axis plotting range
    :param O32_xra: array. x-axis plotting range
    :param yra: array. y-axis plotting range
    :param ctype: list. list of colors
    :param marker: list. list of markers
    :param label: list. list of string labels
    :param log: LogClass or logging object

    No returns

    Notes
    -----
    Created by Chun Ly, 29 November 2018
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    fig, ax = plt.subplots(ncols=2)

    n_sample = len(lR23)

    # Creating the inning by lO32 values for plots
    bin_start = np.zeros(n_bins)
    bin_end = np.zeros(n_bins)

    sort0 = np.argsort(lO32)
    y_sort0 = lO32[sort0]
    r_bin_pts = np.int(np.round(len(lO32) / float(n_bins)))
    lo32_bin_avg = np.zeros(n_bins)

    min1, max1 = np.zeros(n_sample), np.zeros(n_sample)
    OH_min1, OH_max1 = np.zeros(n_sample), np.zeros(n_sample)
    for nn in range(n_sample):
        min1[nn] = np.min(lO32[nn])
        max1[nn] = np.max(lO32[nn])

        OH_min1[nn] = np.min(OH[nn])
        OH_max1[nn] = np.max(OH[nn])

    O32_min = np.min(min1)
    O32_max = np.max(max1)
    O32_arr = np.arange(O32_min, O32_max, 0.025)

    bian_OH = bian18_OH_O32(O32_arr)

    # Grid of 12+log(O/H)
    if len(yra) == 0:
        log.info(len(yra))
        log.info(min(OH_min1))
        OH_max = np.nanmax(OH_max1)
        log.info(OH_max)
        OH_arr = np.arange(min(OH_min1), OH_max, 0.25)
    else:
        OH_arr = np.arange(yra[0], yra[1], 0.05)
    log.info(f"OH_arr: {OH_arr}")
    bian_R23 = bian18_R23_OH(OH_arr, bian_coeff)

    if len(ctype) == 0:
        ctype = ['blue'] * n_sample
    if len(ctype) == 0:
        label = [''] * n_sample

    for nn in range(n_sample):
        for ii in range(n_bins):
            if ii == 0:
                bin_start[0] = y_sort0[0]
                bin_end[0] = y_sort0[r_bin_pts - 1]
            else:
                bin_start[ii] = bin_end[ii - 1] + 0.000001
                bin_end[ii] = y_sort0[
                    np.min([len(lO32) - 1, (ii + 1) * r_bin_pts - 1])]
            print('bin start: ', bin_start[ii], 'bin end: ', bin_end[ii])
            lo32_in_bin = np.where((y_sort0 >= bin_start[ii]) &
                                   (y_sort0 < bin_end[ii]))[0]
            # Plotting
            ax[0].scatter(lR23[nn], OH[nn], color=ctype[ii], edgecolor='none',
                          marker=marker[nn], alpha=0.5, label=label[nn])
            if len(ID[nn]) != 0:
                for aa in range(len(lR23[nn])):
                    ax[0].annotate(ID[nn][aa], [lR23[nn][aa], OH[nn][aa]],
                                   xycoords='data', size=4, va='bottom', ha='left')

            if len(OH_err) != 0:
                ax[0].errorbar(lR23[nn], OH[nn], yerr=OH_err[nn],mec=ctype[ii],
                               ecolor=ctype[nn], capsize=0, alpha=0.5,
                               fmt=None, label=None)

        if len(lR23_err) != 0:
            ax[0].errorbar(lR23[nn], OH[nn], xerr=lR23_err[nn], mec=ctype[ii],
                           ecolor=ctype[nn], capsize=0, alpha=0.5, fmt=None,
                           label=None)

    ax[0].plot(bian_R23, OH_arr, 'k--')

    avail_idx = np.where((OH_arr >= 7.80) & (OH_arr <= 8.4))[0]
    ax[0].plot(bian_R23[avail_idx], OH_arr[avail_idx], 'k-',
               label='Bian+(2018)')

    if len(R23_xra) != 0:
        ax[0].set_xlim(R23_xra)
    if len(yra) != 0:
        ax[0].set_ylim(yra)

    ax[0].set_xlabel(r'$\log(R_{23})$')
    ax[0].set_ylabel(r'$12+\log({\rm O/H})_{T_e}$')
    ax[0].legend(loc='lower left', framealpha=0.5, fontsize=10)

    for nn in range(n_sample):
        ax[1].scatter(lO32[nn], OH[nn], color=ctype[nn], edgecolor='none',
                      marker=marker[nn], alpha=0.5)
        if len(ID[nn]) != 0:
            for ii in range(len(lO32[nn])):
                ax[1].annotate(ID[nn][ii], [lO32[nn][ii], OH[nn][ii]],
                               xycoords='data', size=4, va='bottom', ha='left')

        if len(OH_err) != 0:
            ax[1].errorbar(lO32[nn], OH[nn], yerr=OH_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5, fmt=None,
                           label=None)

        if len(lO32_err) != 0:
            ax[1].errorbar(lO32[nn], OH[nn], xerr=lO32_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5, fmt=None,
                           label=None)

    ax[1].plot(O32_arr, bian_OH, 'k--', label='Bian+(2018)')

    if len(O32_xra) != 0:
        ax[1].set_xlim(O32_xra)
    if len(yra) != 0:
        ax[1].set_ylim(yra)

    ax[1].set_yticklabels([])
    ax[1].set_xlabel(r'$\log(O_{32})$')

    plt.subplots_adjust(left=0.065, right=0.99, bottom=0.1, top=0.99,
                        wspace=0.01)
    fig.set_size_inches(10, 5)
    fig.savefig(out_pdf, format='pdf')

    plot_difference_curvefit.\
        plot_difference_twovariable(lR23, lO32, OH, bin_start, bin_end,
                                    R23_pdf_file, data_input='R23',
                                    new_coefficients=[], n_bins=4,
                                    data_err=lR23_err, OH_err=OH_err,
                                    OH_range=yra, data_range=[-0.5, 0.5],
                                    marker=marker, label=label, IDs=ID, log=log)

    log.info("finished.")


def get_measurements(data, log=None):
    """
    Used in get_DEEP2 and get_MACT to pull data and return it to
    DEEP2_OIII4363 and MACT_OIII4363
    Used in main()
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    lR23 = np.log10(data['R23'].data)
    lO32 = np.log10(data['O32'].data)
    OH = data['OH'].data

    OH_err = np.row_stack((data['OH_lo'].data, data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo, lR23_hi))

    lO32_lo = lO32 - np.log10(data['O32'].data - data['O32_lo'].data)
    lO32_hi = np.log10(data['O32'].data + data['O32_hi'].data) - lO32
    lO32_err = np.row_stack((lO32_lo, lO32_hi))

    log.info("finished.")
    return lR23, lO32, OH, OH_err, lR23_err, lO32_err


def get_DEEP2(path0, log=None):
    """
    Called by DEEP2_OIII4363 and returns data to main function
    Used in main()
    """
    infile = path0 + 'DEEP2_R23_O32_derived.tbl'

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    log.info(f"Reading: {infile}")
    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_measurements(data)

    log.info("finished.")
    return data, lR23, lO32, OH, OH_err, lR23_err, lO32_err


def get_MACT(path0, log=None):
    """
    Called by MACT_OIII4363 and returns data to main function
    Used in main()
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    infile = path0 + 'MACT_R23_O32_derived.tbl'

    log.info(f"Reading : {infile}")

    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_measurements(data)

    log.info("finished.")
    return data, lR23, lO32, OH, OH_err, lR23_err, lO32_err


def DEEP2_OIII4363():
    """
    Run function for DEEP2 dataset for hardcoded path (line 241)
    """

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    log = LogClass(path0,
                   'local_analog_calibration_deep2_oiii4363.log').get_logger()

    log.info("starting ...")

    data, lR23, lO32, OH, OH_err, lR23_err, \
        lO32_err = get_DEEP2(path0, log=log)

    ID = data['ID']

    out_pdf = path0 + 'DEEP2_R23_O32_Bian18.pdf'
    main([lR23], [lO32], [OH], out_pdf, ID=[ID], lR23_err=[lR23_err],
         lO32_err=[lO32_err], OH_err=[OH_err], R23_xra=[0.75, 1.05],
         O32_xra=[0.05, 0.95], yra=[7.1, 8.65], label=[r'DEEP2 (Ly+2015)'],
         log=log)

    log.info("finished.")


def MACT_OIII4363():
    """
    Run function for MACT dataset for hardcoded path (line 257) 
    """
    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    log = LogClass(path0,
                   'local_analog_calibration_mact_oiii4363.log').get_logger()

    log.info("starting ...")

    data, lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_MACT(path0, log=log)

    ID = data['ID']

    out_pdf = path0 + 'MACT_R23_O32_Bian18.pdf'
    main([lR23], [lO32], [OH], out_pdf, ID=[ID], lR23_err=[lR23_err],
         lO32_err=[lO32_err], OH_err=[OH_err], R23_xra=[0.6, 1.15],
         O32_xra=[-0.55, 2.1], yra=[7.1, 8.85],
         label=[r'$\mathcal{MACT}$  (Ly+2016)'], log=log)

    log.info("finished.")


def DEEP2_MACT_OIII4363():
    """
    Run function for DEEP2 and MACT (combined) dataset for hardcoded path
    (line 277)
    """
    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    log = LogClass(path0,
                   'local_analog_calibration_deep2_mact_oiii4363.log').get_logger()

    log.info("starting ...")

    # DEEP2
    DEEP2_data, DEEP2_lR23, DEEP2_lO32, DEEP2_OH, DEEP2_OH_err, \
        DEEP2_lR23_err, DEEP2_lO32_err = get_DEEP2(path0, log=log)

    # MACT
    MACT_data, MACT_lR23, MACT_lO32, MACT_OH, MACT_OH_err, \
        MACT_lR23_err, MACT_lO32_err = get_MACT(path0, log=log)

    # Combine together
    lR23 = [DEEP2_lR23, MACT_lR23]
    lO32 = [DEEP2_lO32, MACT_lO32]
    OH = [DEEP2_OH,   MACT_OH]

    OH_err = [DEEP2_OH_err, MACT_OH_err]

    lR23_err = [DEEP2_lR23_err, MACT_lR23_err]
    lO32_err = [DEEP2_lO32_err, MACT_lO32_err]

    ID = [DEEP2_data['ID'].data, MACT_data['ID'].data]

    out_pdf = path0 + 'DEEP2_MACT_R23_O32_Bian18.pdf'

    labels = [r'DEEP2 [OIII]$\lambda$4363-detected',
              r'$\mathcal{MACT}$  (Ly+2016)']
    main(lR23, lO32, OH, out_pdf, ID=ID, lR23_err=lR23_err, lO32_err=lO32_err,
         OH_err=OH_err, R23_xra=[0.6, 1.15], O32_xra=[-0.55, 2.1],
         yra=[7.1, 8.85], ctype=['blue', 'green'], label=labels, log=log)

    out_pdf = path0 + 'DEEP2_MACT_R23_O32_Bian18.nolabel.pdf'

    main(lR23, lO32, OH, out_pdf, lR23_err=lR23_err, lO32_err=lO32_err,
         OH_err=OH_err, R23_xra=[0.6, 1.15], O32_xra=[-0.55, 2.1],
         yra=[7.1, 8.85], ctype=['blue', 'green'], label=labels,
         log=log)

    log.info("finished.")


def zcalbase():
    """
    Run function for Zcalbase_gal dataset for hardcoded path (line 316) 
    """
    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    log = LogClass(path0, 'local_analog_calibration_zcalbase.log').get_logger()

    log.info("starting ...")

    ref_name0 = ['Berg2012', 'Kennicutt2003', 'Izotov1994', 'Thuan1995',
                 'Izotov1997', 'Guseva2009', 'Izotov2012', 'SDSS']
    dir0 = ['', '', 'BCGs', 'BCGs', 'Pilyugin2012/Izotov1997',
            'Pilyugin2012/Guseva2009', 'Pilyugin2012/Izotov2012', '']

    lR23_all = []
    lO32_all = []
    OH_all = []

    for name, dir_path in zip(ref_name0, dir0):
        lR23, lO32, OH = get_zcalbase_sample(name, dir_path=dir_path,
                                             log=log)

        lR23_all.append(lR23)
        lO32_all.append(lO32)
        OH_all.append(OH)

    out_pdf = path0 + 'Zcalbase_Bian18.pdf'
    label = ref_name0  # ['Kennicutt+2003']
    ctype = ['blue', 'cyan', 'green', 'yellow',
             'red', 'magenta', 'gray', 'black']
    main(lR23_all, lO32_all, OH_all, out_pdf, R23_xra=[0.1, 1.05],
         O32_xra=[-0.5, 1.5], yra=[7.0, 9.0], label=label,
         ctype=ctype, log=log)
