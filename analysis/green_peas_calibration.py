
from chun_codes import match_nosort_str
from astropy.io import ascii as asc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.io import fits
from os.path import join

from .deep2_r23_o32.log_commons import log_stdout, LogClass

jiang18_coeffs = [-24.135, 6.1523, -0.37866, -0.147, -7.071]


def O32_OH_fit(xy, a, b, c, d, e):
    """
    Main functional code that determine log(R23) from log(O32) and 12+log(O/H)
    Used by main() function
    :param x: 12+log(O/H)
    :param y: log([OIII]/[OII])
    :param a, b, c, d, e,:  Jiang coefficients
    """

    x = xy[0]
    y = xy[1]

    logR23 = a + b * x + c * x**2 - d * (e + x) * y

    return logR23


def jiang18(x, y):
    """
    Function to return log(R23) based on metallicity and [OIII]/[OII] flux ratio
    Used by main() function

    :param x: array. 12+log(O/H)
    :param y: array. log([OIII]/[OII])
    """

    logR23 = O32_OH_fit(xy, * jiang18_coeffs)

    return logR23


def plot_differences(lR23, lO32, OH, lO32_all, out_diff_pdf, bin_start,
                     bin_end, n_bins=4, lR23_err=[], OH_err=[], OH_range=[],
                     dR23_range=[], marker=[], label=[], IDs=[], log=None):
    """
    Plot differences between Jiang18 R23 vs observed R23 as a function of metallicity
    Used by main() function
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    # Plotting
    ctype = ['red', 'magenta', 'green', 'cyan', 'blue', 'black']

    if len(marker) == 0:
        marker = ['o'] * n_sample

    diff0 = []
    for nn in range(n_sample):
        jiang_R23 = O32_OH_fit((OH[nn], lO32[nn]), *jiang18_coeffs)
        log.info(f"jiang_R23: {jiang_R23}")

        if nn == 0:
            if IDs:
                for jj in range(len(jiang_R23)):
                    id_diff = lR23[0][jj]-jiang_R23[jj]
                    ax.annotate(IDs[0][jj], (OH[0][jj], id_diff), fontsize='6')

        # Label in upper left the points
        if len(label) != 0:
            x1 = OH_range[0] + 0.025*(OH_range[1]-OH_range[0])
            y1 = dR23_range[1] - (nn*0.035 + 0.05)*(dR23_range[1]-dR23_range[0])
            x2 = OH_range[0] + 0.035*(OH_range[1]-OH_range[0])
            y2 = dR23_range[1] - (nn*0.035 + 0.0525)*(dR23_range[1]-dR23_range[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1], [y1], marker=marker[nn], color='black')

        for ii in range(n_bins):
            y_ii_min = bin_start[ii]
            y_ii_max = bin_end[ii]
            idx = np.where((lO32[nn] >= y_ii_min) & (lO32[nn] <= y_ii_max))[0]

            ii_label = ''
            if nn == 0:  # n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) & (lO32_all <= y_ii_max))[0]
                ii_label = fr" {y_ii_min:.2f} < $\log(O_{{32}})$ " + \
                           fr"< {y_ii_max:.2f}, N = {len(idx_all):d}"

            if len(idx) > 0:
                i_diff = lR23[nn][idx] - jiang_R23[idx]
                ax.scatter(OH[nn][idx], i_diff, color=ctype[ii], marker=marker[nn],
                           edgecolor='none', alpha=0.5, label=ii_label)

                diff0 += list(lR23[nn][idx] - jiang_R23[idx])

                # Added if statement so that only data points on the OH_err[0] place will be plotted
                if nn == 0: 
                    if len(OH_err) != 0:
                        ax.errorbar(OH[nn][idx], i_diff, xerr=np.transpose(OH_err[nn][idx]),
                                    mfc='none', ecolor=ctype[ii], capsize=0, alpha=0.25,
                                    fmt='None', label=None, ls='none')

                    if len(lR23_err) != 0:
                        ax.errorbar(OH[nn][idx], i_diff, yerr=np.transpose(lR23_err[nn][idx]),
                                    mfc='none', ecolor=ctype[ii], capsize=0, alpha=0.25,
                                    fmt='None', label=None, ls='none')

    # Draw horizontal line at zero:
    ax.axhline(y=0, c='k', linestyle='dashed')

    # Compute statistics
    med0 = np.median(diff0)
    avg0 = np.average(diff0)
    sig0 = np.std(diff0)
    ax.axhline(y=avg0, c='r', linestyle='dotted')
    ax.axhline(y=med0, c='b', linestyle='dotted')

    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.annotate(an_txt, [0.2, 0.015], xycoords='axes fraction', va='bottom', ha='right',
                fontsize=10)

    if len(OH_range) != 0:
        ax.set_xlim(OH_range)
    if len(dR23_range) != 0:
        ax.set_ylim(dR23_range)

    ax.set_xlabel(r'$12+\log({\rm O/H})_{T_e}$')
    ax.set_ylabel(r'$\Delta_{R_{23}} \equiv \log(R_{23}) - \log(R_{23})_{\rm GPC}$')
    ax.minorticks_on()

    leg = ax.legend(loc='upper right', scatterpoints=1, fontsize=8, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.1, top=0.97)

    log.info(f"Writing: {out_diff_pdf}")
    fig.savefig(out_diff_pdf)

    log.info("finished.")


def main(lR23, lO32, OH, out_pdf, n_bins=4, lR23_err=[], OH_err=[], xra=[], yra=[],
         marker=[], edgecolors=[], alpha=[], label=[], dR23_range=[-0.3, 0.3],
         IDs=[], include_Rlimit=True, fit=False, log=None):
    """
    Plots R23 and O32 line ratios and compare it to Jiang et al. (2018)
    calibration that is based on green pea galaxies
    Main function to plot dataset against Jiang+ (2018) calibration

    :param lR23: array. log(R_23)
    :param lO32: array. log([OIII]/[OII])
    :param OH: array. 12+log(O/H)
    :param out_pdf: str. full path for output PDF
    :param n_bins: int. number of bins to bin data
    :param lR23_err: array. error on x-axis (logR23)
    :param OH_err: array. error on y-axis (logR23)
    :param xra: array. x-axis plotting range
    :param yra: array. y-axis plotting range
    :param marker: list. list of markers
    :param edgecolors: list. list of edge colors
    :param alpha: list. list of alpha (transparency)
    :param label: list. list of string labels
    :param dR23_range: array. for annotation
    :param IDs: array. list of list of ID's
    :param include_Rlimit: Bool. to plot cases with R-limits
    :param fit: bool. Turn on fitting. Default: False -> Uses Jiang+2018 relation
    :param log: LogClass or logging object

    No returns

    Notes
    -----
    Created by Chun Ly, 23 November 2018

    Edited by Reagen Leimbach, 18 June 2020 to add IDs to the plots
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    min1, max1 = np.zeros(n_sample), np.zeros(n_sample)
    OH_min1, OH_max1 = np.zeros(n_sample), np.zeros(n_sample)
    for nn in range(n_sample):
        min1[nn] = np.min(lO32[nn])
        max1[nn] = np.max(lO32[nn])

        good = np.where(np.isfinite(OH[nn]))[0]
        OH_min1[nn] = np.min(OH[nn][good])
        OH_max1[nn] = np.max(OH[nn][good])

    # Grid of 12+log(O/H)
    if len(yra) == 0:
        x_arr = np.arange(min(OH_min1), max(OH_max1), 0.05)
    else:
        x_arr = np.arange(yra[0], yra[1], 0.05)

    lO32_all = np.array([])
    for nn in range(n_sample):
        lO32_all = np.append(lO32_all, lO32[nn])

    # Binning based on O32 data because of the way the Jiang calibration works
    sort0 = np.argsort(lO32_all)
    y_sort0 = lO32_all[sort0]

    # bin_pts = np.int(len(lO32_all)/n_bins)
    r_bin_pts = np.int(np.round(len(lO32_all)/float(n_bins)))

    bin_start = np.zeros(n_bins)
    bin_end = np.zeros(n_bins)
    
    bin_start[0] = y_sort0[0]
    bin_end[0] = y_sort0[r_bin_pts-1]
    for ii in range(1, n_bins):
        bin_start[ii] = bin_end[ii-1]+0.000001
        bin_end[ii] = y_sort0[np.min([len(lO32_all)-1, (ii+1)*r_bin_pts-1])]

    # Plotting
    ctype = ['red', 'magenta', 'green', 'cyan', 'blue', 'black']
    xytext_location = ([5, 2], [5, 2], [0, 2], [3, 2], [7, 2], [5, 2])

    if len(marker) == 0:
        marker = ['o'] * n_sample

    if fit:
        p0 = jiang18_coeffs

        OH_all = np.array([])
        lR23_all = np.array([])
        for nn in range(n_sample):
            OH_all = np.append(OH_all, OH[nn])
            lR23_all = np.append(lR23_all, lR23[nn])

        opt, cov = curve_fit(O32_OH_fit, (OH_all, lO32_all), lR23_all, p0=p0)
        log.info(opt)

    for nn in range(n_sample):
        if len(label) != 0:
            x1 = xra[0] + 0.28*(xra[1]-xra[0])
            y1 = yra[1] - (nn*0.035 + 0.1)*(yra[1]-yra[0])
            x2 = xra[0] + 0.3*(xra[1]-xra[0])
            y2 = yra[1] - (nn*0.035 + 0.1)*(yra[1]-yra[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1], [y1], marker=marker[nn], color='black')

        for ii in range(n_bins):
            y_ii_min = bin_start[ii]  # bin_y_min + ii * dy
            y_ii_max = bin_end[ii]    # y_min + (ii+1) * dy
            idx = np.where((lO32[nn] >= y_ii_min) & (lO32[nn] <= y_ii_max))[0]
            log.info(f"idx: {idx}")

            ii_label = ''
            if nn == 0:  # n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) & (lO32_all <= y_ii_max))[0]
                ii_label = fr" {y_ii_min:.2f} < $\log(O_{{32}})$ " + \
                           f"< {y_ii_max:.2f}, N = {len(idx_all):d}"
            if len(idx) > 0:
                ax.scatter(lR23[nn][idx], OH[nn][idx], color=ctype[ii], marker=marker[nn],
                           alpha=alpha[nn], label=ii_label, edgecolors=edgecolors[nn])

                # Pushed Error bars under idx requirement
                # Added if statement so that only data points on the OH_err[0] place will be plotted
                if nn == 0:
                    if len(OH_err) != 0:
                        log.info(f"OH: {OH[nn][idx]}")
                        log.info(f"OH_err: {OH_err[nn][idx]}")
                        ax.errorbar(lR23[nn][idx], OH[nn][idx], yerr=np.transpose(OH_err[nn][idx]),
                                    mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5,
                                    fmt='', label=None, ls='none')

                    if len(lR23_err) != 0:
                        ax.errorbar(lR23[nn][idx], OH[nn][idx], xerr=np.transpose(lR23_err[nn][idx]),
                                    mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5, fmt=None,
                                    label=None, ls='none')

            if nn == 0:  # n_sample-1:
                lO32_avg = np.average(lO32_all[idx_all])
                if not fit:
                    opt = jiang18_coeffs

                mod_logR23 = O32_OH_fit((x_arr, lO32_avg), *opt)
                ax.annotate(f"{lO32_avg:.2f}", [mod_logR23[-1], x_arr[-1]],
                            xytext=xytext_location[ii], textcoords='offset points',
                            color=ctype[ii], xycoords='data', ha='center',
                            va='bottom', fontsize=8)
                ax.plot(mod_logR23, x_arr, color=ctype[ii], linestyle='dashed')

    ax.set_xlim(0.5, 1.1)
    ax.set_ylim(6.75, np.max(x_arr) + 0.1)
    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})_{T_e}$')
    leg = ax.legend(loc='lower right', scatterpoints=1, fontsize=8,
                    framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

    # This puts the IDs on all the given IDs entered into the R23 and OH arrays
    if IDs: 
        for yy in range(len(IDs)):
            id_a = IDs[yy]
            R23_a = lR23[yy]
            OH_a = OH[yy]
            for aa in range(len(id_a)):
                ax.annotate(id_a[aa], (R23_a[aa], OH_a[aa]), fontsize='6')
                
    plt.subplots_adjust(left=0.075, right=0.99, bottom=0.08, top=0.97)
    fig.savefig(out_pdf)

    # Because we do not want to include the Robust limits into the statistical calculations
    # in plot_differences, this options allows to redefine lR23, lO32, OH
    if include_Rlimit:
        nR23 = [lR23[0], lR23[2], lR23[3]]
        nO32 = [lO32[0], lO32[2], lO32[3]]
        nOH = [OH[0], OH[2], OH[3]]
        nIDs = [IDs[0]]
        log.info('Using redefined values')
        label = ['Detection', 'DEEP2', 'MACT']
        marker = ['D', '3', '4']
    else:
        nR23 = lR23
        nO32 = lO32
        nOH = OH
        nIDs = IDs

    # Plot differences between model and data
    if not fit:
        out_diff_pdf = out_pdf.replace('.pdf', '.diff.pdf')
        plot_differences(nR23, nO32, nOH, lO32_all, out_diff_pdf, bin_start, bin_end,
                         n_bins=n_bins, lR23_err=lR23_err, OH_err=OH_err,
                         OH_range=yra, dR23_range=dR23_range, marker=marker,
                         label=label, IDs=nIDs, log=log)

    log.info("finished.")


def get_measurements(data, log=None):
    """
    Used in get_DEEP2 and get_MACT to pull data and return it to
    DEEP2_OIII4363 and MACT_OIII4363
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    lR23 = np.log10(data['R23'].data)
    lO32 = np.log10(data['O32'].data)
    OH = data['OH'].data

    OH_err = np.row_stack((data['OH_lo'].data, data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo, lR23_hi))

    log.debug("finished.")

    return lR23, lO32, OH, OH_err, lR23_err


def get_zcalbase_sample(prefix, dir_path='', log=None):
    """
    Used in Zcalbase function to pull data
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    dir0 = '/Users/cly/data/Metallicity/Others/Te_Repository/'
    if dir_path == '':
        flux_file = f"{dir0}{prefix}/{prefix}_sample.det4363.int.fits"
        Te_file = f"{dir0}{prefix}/{prefix}_Te_table.fits"
    else:
        flux_file = f"{dir0}{dir_path}/{prefix}_sample.det4363.int.fits"
        Te_file = f"{dir0}{dir0}/{prefix}_Te_table.fits"

    if 'SDSS' in prefix:
        path_SDSS = '/Users/cly/data/Metallicity/Others/SDSS/'
        flux_file = join(path_SDSS, 'SDSS_DR7_det4363.OIIfix.int.fits')
        Te_file = join(path_SDSS, 'SDSS_DR7_det4363_Te_table.dered.fits')

    data0 = fits.getdata(flux_file)
    data1 = fits.getdata(Te_file)

    lR23 = np.log10((data0.OII_3727_FLUX + data0.OIII_5007_FLUX) / data0.H_BETA_FLUX)
    lO32 = np.log10(data0.OIII_5007_FLUX / data0.OII_3727_FLUX)

    OH = data1.OH_gas

    log.info("finished.")

    return lR23, lO32, OH


def get_DEEP2(path0, log=None):
    """
    Called by DEEP2_OIII4363 and returns data to main function
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    infile = join(path0, 'DEEP2_R23_O32_derived.tbl')

    log.info(f"Reading: {infile}")
    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err = get_measurements(data, log=log)

    log.debug("finished.")

    return data, lR23, lO32, OH, OH_err, lR23_err


def get_MACT(path0, log=None):
    """
    Called by MACT_OIII4363 and returns data to main function
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    infile = join(path0, 'MACT_R23_O32_derived.tbl')

    log.info(f"Reading: {infile}")
    data = asc.read(infile)

    dup = ['Keck10', 'Keck17', 'Keck22', 'Keck25']
    d_idx1, d_idx2 = match_nosort_str(data['ID'].data, dup)

    data.remove_rows(d_idx1)

    lR23, lO32, OH, OH_err, lR23_err = get_measurements(data, log=log)

    log.debug("finished.")

    return data, lR23, lO32, OH, OH_err, lR23_err


def DEEP2_OIII4363(log_dir):
    """
    Run function for DEEP2 dataset for hardcoded path (line 442)
    """

    log = LogClass(log_dir, 'green_peas_calibration_deep2_oiii4363.log').get_logger()

    log.debug("starting ...")

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err = get_DEEP2(path0, log=log)

    out_pdf = join(path0, 'DEEP2_R23_O32_Jiang18.pdf')
    main([lR23], [lO32], [OH], out_pdf, lR23_err=[lR23_err], OH_err=[OH_err],
         xra=[0.75, 1.05], yra=[7.1, 8.65], log=log)

    # Got RuntimeError
    # out_pdf = join(path0, 'DEEP2_R23_O32_Jiang18.fit.pdf')
    # main([lR23], [lO32], [OH], out_pdf, lR23_err=[lR23_err], OH_err=[OH_err],
    # xra=[0.75,1.05], yra=[7.1,8.65], fit=True)

    log.debug("finished.")


def MACT_OIII4363(log_dir):
    """
    Run function for MACT dataset for hardcoded path (line 462)
    """

    log = LogClass(log_dir, 'green_peas_calibration_mact_oiii4363.log').get_logger()

    log.debug("starting ...")

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err = get_MACT(path0, log=log)

    out_pdf = join(path0, 'MACT_R23_O32_Jiang18.pdf')
    main([lR23], [lO32], [OH], out_pdf, n_bins=6, lR23_err=[lR23_err],
         OH_err=[OH_err], xra=[0.60, 1.15], yra=[7.10, 8.7], log=log)

    out_pdf = join(path0, 'MACT_R23_O32_Jiang18.fit.pdf')
    main([lR23], [lO32], [OH], out_pdf, n_bins=6, lR23_err=[lR23_err],
         OH_err=[OH_err], xra=[0.60, 1.15], yra=[7.10, 8.7], fit=True,
         log=log)

    log.debug("finished.")


def DEEP2_MACT_OIII4363(log_dir, include_stack=False, fit=False):
    """
    Run function for DEEP2 and MACT (combined) dataset for hardcoded path (line 481)
    """

    log = LogClass(log_dir, 'green_peas_calibration_deep2_mact_oiii4363.log').get_logger()

    log.debug("starting ...")

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    # DEEP2
    DEEP2_data, DEEP2_lR23, DEEP2_lO32, DEEP2_OH, DEEP2_OH_err, \
        DEEP2_lR23_err = get_DEEP2(path0, log=log)

    # MACT
    MACT_data, MACT_lR23, MACT_lO32, MACT_OH, MACT_OH_err, \
        MACT_lR23_err = get_MACT(path0, log=log)

    # Combine together
    lR23 = [MACT_lR23, DEEP2_lR23]
    lO32 = [MACT_lO32, DEEP2_lO32]
    OH = [MACT_OH,   DEEP2_OH]

    OH_err = [MACT_OH_err, DEEP2_OH_err]
    lR23_err = [MACT_lR23_err, DEEP2_lR23_err]

    # Include RLeimbach stacked results
    if include_stack:
        stack_infile = '/Users/cly/Google Drive/Zcalbase_gal/' + \
                        'Double_Bin_temperatures_metalicity.tbl'
        log.info('### Reading : ' + stack_infile)
        stack_tbl = asc.read(stack_infile)

        stack_det = np.where((stack_tbl['S/N_4363'] >= 3.0) &
                             (stack_tbl['com_O_log'] >= 8.1))[0]
        log.info(f"stack_det: {len(stack_det)}")
        R23_stack = stack_tbl['R23_Composite'].data[stack_det]
        O32_stack = stack_tbl['O32_Composite'].data[stack_det]
        OH_stack = stack_tbl['com_O_log'].data[stack_det]

        lR23 += [R23_stack]
        lO32 += [O32_stack]
        OH += [OH_stack]

        OH_err += [np.array([[0] * len(stack_det), [0] * len(stack_det)])]
        lR23_err += [np.array([[0] * len(stack_det), [0] * len(stack_det)])]
        log.info(OH_err)

    out_pdf = join(path0, 'MACT_DEEP2_R23_O32_Jiang18.pdf')
    label = [r'$\mathcal{MACT}$  (Ly+2016)',
             r'DEEP2 [OIII]$\lambda$4363-detected']
    marker = ['o', '*']

    if include_stack:
        out_pdf = out_pdf.replace('.pdf', '.stack.pdf')
        label += ['DEEP2 stacked detections']
        marker += ['s']
        yra = [7.10, 9.0]
    else:
        yra = [7.10, 8.8]

    main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.45, 1.15], yra=yra, dR23_range=[-0.15, 0.3], marker=marker,
         label=label, log=log)

    if fit:
        out_pdf = join(path0, 'MACT_DEEP2_R23_O32_Jiang18.fit.pdf')
        main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
             xra=[0.45, 1.15], yra=yra, marker=['*', 'o'], label=label,
             fit=True, log=log)

    log.debug("finished.")


def zcalbase(log_dir):
    """
    Run function for Zcalbase_gal dataset for hardcoded path (line 548)
    """

    log = LogClass(log_dir, 'green_peas_calibration_zcalbase.log').get_logger()

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    ref_name0 = ['Berg2012', 'Kennicutt2003', 'Izotov1994', 'Thuan1995',
                 'Izotov1997', 'Guseva2009', 'Izotov2012', 'SDSS']
    dir0 = ['', '', 'BCGs', 'BCGs', 'Pilyugin2012/Izotov1997',
                'Pilyugin2012/Guseva2009', 'Pilyugin2012/Izotov2012', '']

    lR23_all = []
    lO32_all = []
    OH_all = []

    for name, dir0 in zip(ref_name0, dir0):
        lR23, lO32, OH = get_zcalbase_sample(name, dir_path=dir0, log=log)

        lR23_all.append(lR23)
        lO32_all.append(lO32)
        OH_all.append(OH)

    out_pdf = join(path0, 'Zcalbase_Jiang18.pdf')
    label = ref_name0  # ['Kennicutt+2003']
    main(lR23_all, lO32_all, OH_all, out_pdf, n_bins=6, xra=[0.1, 1.05],
         yra=[7.0, 9.0], marker=['s', 'x', 's', 's', 'D', 'D', 's', '*'],
         label=label, log=log)
