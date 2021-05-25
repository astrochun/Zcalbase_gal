import numpy as np
from astropy.io import ascii as asc
from scipy.optimize import curve_fit
from os.path import join
import matplotlib.pyplot as plt

from Zcalbase_gal.analysis.deep2_r23_o32.log_commons \
    import log_stdout, LogClass
from Zcalbase_gal.analysis.deep2_r23_o32 import bian_coeff, ctype, \
    secondorder_polynomial, thirdorder_polynomial, threevariable_fit, \
    bian18_R23_OH, bian18_O32_OH, bian18_OH_O32, jiang_O32_OH_fit, jiang18_coeffs


def plot_difference_threevariable(lR23, lO32, OH, lO32_all, pdf_file, fitting_model,
                                  bin_start, bin_end, new_coefficients=[],
                                  n_bins=4, data_err=[], OH_err=[],
                                  OH_range=[7,11], data_range=[-0.3, 0.3],
                                  marker=[], label=[], IDs=[], log=None):
    """
    So this curve fitting function looks like it fits for three variables fits
    a * x**3 + b*x**2 + c * x + d*lO32
    This is also going to be used for the Jiang calibration


    Plot differences between LACR23 vs observed R23 as
    a function of metallicity
    Plot differences between LACO32 vs observed O32 as
    a function of metallicity
    Used by main() function
    """
    # print('new coeff: ', new_coefficients, type(new_coefficients))
    if log is None:
        log = log_stdout()

    # log.info("starting ...")

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    # Plotting
    if len(marker) == 0:
        marker = ['o'] * n_sample

    diff0 = []
    for nn in range(n_sample):
        if fitting_model == 'jiang':
            fitted_function = jiang_O32_OH_fit((OH[nn], lO32[nn]),
                                               *jiang18_coeffs)
            log.info(f"jiang_R23: {fitted_function}")
        if fitting_model == 'Zcal':
            if len(new_coefficients) != 0:
                fitted_function = threevariable_fit((OH[nn], lO32[nn]),
                                                    *new_coefficients)
                # log.info(f"curve fit LAC_R23: {fitted_function}")

        if nn == 0:
            if IDs:
                for jj in range(len(fitted_function)):
                    id_diff = lR23[nn][jj] - fitted_function[jj]
                    ax.annotate(IDs[nn][jj], (OH[nn][jj], id_diff),
                                fontsize='6')

        # Label in upper left the points
        if len(label) != 0:
            x1 = OH_range[0] + 0.025 * (OH_range[1] - OH_range[0])
            y1 = data_range[1] - (nn * 0.035 + 0.05) \
                     * (data_range[1] - data_range[0])
            x2 = OH_range[0] + 0.035 * (OH_range[1] - OH_range[0])
            y2 = data_range[1] - (nn * 0.035 + 0.0525) \
                     * (data_range[1] - data_range[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1], [y1], marker=marker[nn], color='black')

        for ii in range(n_bins):
            y_ii_min = bin_start[ii]
            y_ii_max = bin_end[ii]
            idx = np.where((lO32[nn] >= y_ii_min) &
                           (lO32[nn] <= y_ii_max))[0]

            ii_label = ''
            if nn == 0:  # n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) &
                                   (lO32_all <= y_ii_max))[0]
                ii_label = fr" {y_ii_min:.2f} < $\log(O_{{32}})$ " + \
                           fr"< {y_ii_max:.2f}, N = {len(idx_all):d}"

            if len(idx) > 0:
                i_diff = lR23[nn][idx] - fitted_function[idx]
                ax.scatter(OH[nn][idx], i_diff, color=ctype[nn], marker=marker[nn],
                           edgecolor='none', alpha=0.5, label=ii_label)
                if nn == 1:
                    ax.scatter(OH[nn][idx], i_diff, color=ctype[1],
                               marker=r'$\rightarrow$',
                               edgecolor='none', alpha=0.5)
                diff0 += list(lR23[nn][idx] - fitted_function[idx])

        # Added if statement so that only data points
        # on the OH_err[0] place will be plotted
                if nn == 0:
                    if len(OH_err) != 0:
                        ax.errorbar(OH[nn][idx], i_diff,
                                    xerr=np.transpose(OH_err[nn]),
                                    mfc='none', capsize=0,
                                    alpha=0.25, fmt='None', label=None,
                                    ls='none')

                    if len(data_err) != 0:
                        ax.errorbar(OH[nn][idx], i_diff,
                                    yerr=np.transpose(data_err[nn]),
                                    mfc='none', capsize=0,
                                    alpha=0.25, fmt='None', label=None,
                                    ls='none')

    # Draw horizontal line at zero:
    ax.axhline(y=0, c='k', linestyle='dashed')

    # Compute statistics for R23
    med0 = np.median(diff0)
    avg0 = np.average(diff0)
    sig0 = np.std(diff0)

    # Plotting for R23
    ax.axhline(y=avg0, c='r', linestyle='dotted')
    ax.axhline(y=med0, c='b', linestyle='dotted')

    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.set_ylabel(r'$\Delta_{R_{23}} \equiv \log(R_{23}) '
                  r'- \log(R_{23})_{\rm Zcal_curvefit}$')

    ax.annotate(an_txt, [0.2, 0.015], xycoords='axes fraction',
                va='bottom', ha='right', fontsize=10)
    ax.set_xlabel(r'$12+\log({\rm O/H})_{T_e}$')

    ax.minorticks_on()

    if len(OH_range) != 0:
        ax.set_xlim(OH_range)
    if len(data_range) != 0:
        ax.set_ylim(data_range)
    ax.set_xlim(7, 9.5)
    leg_R23 = ax.legend(loc='upper right', scatterpoints=1, fontsize=8,
                            framealpha=0.5)
    for lh in leg_R23.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.1, top=0.97)

    # log.info(f"Writing: {pdf_file}")
    fig.savefig(pdf_file)

    # log.info("finished.")


def plot_difference_twovariable(lR23_lO32, OH, pdf_file, data_input,
                                 new_coefficients=[], data_err=[], OH_err=[],
                                 OH_range=[], data_range=[-0.3, 0.3], marker=[],
                                 label=[], IDs=[], log=None):
    """
    Plot differences between LACR23 vs observed R23 as a function
    of metallicity
    Plot differences between LACO32 vs observed O32 as a function
    of metallicity
    Used by main() function
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    fig, ax = plt.subplots()

    n_sample = len(lR23_lO32)

    # Plotting
    if len(marker) == 0:
        marker = ['o'] * n_sample

    diff0 = []
    for nn in range(n_sample):
        if len(new_coefficients) != 0:
            # This needs some more work
            if data_input == 'R23':
                LAC = bian18_R23_OH(OH[nn], new_coefficients)
                log.info(f"curve fit LAC_R23: {LAC}")
            else:
                LAC = bian18_O32_OH(OH[nn])
                log.info(f"curve fit LAC_R23: {LAC}")
        else:
            if data_input == 'R23':
                LAC = bian18_R23_OH(OH[nn], bian_coeff)
                log.info(f"LAC_R23: {LAC}")
            else:
                LAC = bian18_O32_OH(OH[nn])
                log.info(f"LAC_R23: {LAC}")

        if nn == 0:
            if IDs:
                for jj in range(len(LAC)):
                    id_diff = lR23_lO32[0][jj]-LAC[jj]
                    ax.annotate(IDs[0][jj], (OH[0][jj], id_diff),
                                fontsize='6')

        # Label in upper left the points
        if len(label) != 0:
            x1 = OH_range[0] + 0.025 * (OH_range[1] - OH_range[0])
            y1 = data_range[1] - (nn * 0.035 + 0.05) \
                     * (data_range[1] - data_range[0])
            x2 = OH_range[0] + 0.035 * (OH_range[1] - OH_range[0])
            y2 = data_range[1] - (nn * 0.035 + 0.0525) \
                     * (data_range[1] - data_range[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1], [y1], marker=marker[nn], color='black')
        for ii in range(len(lR23_lO32)):
            y_ii_min = np.min(lR23_lO32[nn])
            y_ii_max = np.max(lR23_lO32[nn])
            idx = np.where((lR23_lO32[nn] >= y_ii_min) &
                           (lR23_lO32[nn] <= y_ii_max))[0]

            if len(idx) > 0:
                i_diff = lR23_lO32[nn][idx] - LAC[idx]
                ax.scatter(OH[nn], i_diff, color=ctype[nn],
                           marker=marker[nn], edgecolor='none', alpha=0.5)
                diff0 += list(lR23_lO32[nn][idx] - LAC[idx])

        # Added if statement so that only data points
        # on the OH_err[0] place will be plotted
                if nn == 0:
                    if len(OH_err) != 0:
                        ax.errorbar(OH[nn], i_diff,
                                    xerr=np.transpose(OH_err[nn]),
                                    mfc='none', capsize=0,
                                    alpha=0.25, fmt='None', label=None,
                                    ls='none')

                    if len(data_err) != 0:
                        ax.errorbar(OH[nn], i_diff,
                                    yerr=np.transpose(data_err[nn]),
                                    mfc='none', capsize=0,
                                    alpha=0.25, fmt='None', label=None,
                                    ls='none')

    # Draw horizontal line at zero:
    ax.axhline(y=0, c='k', linestyle='dashed')

    # Compute statistics for R23
    med0 = np.median(diff0)
    avg0 = np.average(diff0)
    sig0 = np.std(diff0)

    # Plotting for R23
    ax.axhline(y=avg0, c='r', linestyle='dotted')
    ax.axhline(y=med0, c='b', linestyle='dotted')

    if data_input == 'R23':
        an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
        an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
        an_txt += r'$\sigma$ : %0.2f' % sig0
        ax.set_ylabel(r'$\Delta_{R_{23}} \equiv \log(R_{23}) '
                      r'- \log(R_{23})_{\rm LAC}$')
    else:
        an_txt = r'$<\Delta_{O_{32}}>$ : %0.2f' % avg0 + '\n'
        an_txt += r'$\tilde\Delta_{O_{32}}$ : %0.2f' % med0 + '\n'
        an_txt += r'$\sigma$ : %0.2f' % sig0
        ax.set_ylabel(r'$\Delta_{O_{32}} \equiv \log(O_{32}) '
                      r'- \log(O_{32})_{\rm LAC}$')
    ax.annotate(an_txt, [0.2, 0.015], xycoords='axes fraction',
                va='bottom', ha='right', fontsize=10)
    ax.set_xlabel(r'$12+\log({\rm O/H})_{T_e}$')

    ax.minorticks_on()

    if len(OH_range) != 0:
        ax.set_xlim(OH_range)
    if len(data_range) != 0:
        ax.set_ylim(data_range)

    leg_R23 = ax.legend(loc='upper right', scatterpoints=1, fontsize=8,
                            framealpha=0.5)
    for lh in leg_R23.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.1, top=0.97)

    log.info(f"Writing: {pdf_file}")
    fig.savefig(pdf_file)

    log.info("finished.")

