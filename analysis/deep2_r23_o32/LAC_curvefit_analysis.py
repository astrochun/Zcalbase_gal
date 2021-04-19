import numpy as np
from astropy.io import ascii as asc
from scipy.optimize import curve_fit
from os.path import join
import matplotlib.pyplot as plt
from Zcalbase_gal.analysis import local_analog_calibration
from Metallicity_Stack_Commons.column_names import filename_dict
from Zcalbase_gal.analysis.deep2_r23_o32 import bian_coeff
from .log_commons import log_stdout, LogClass


def secondorder_polynomial(x, a, b, c):
    return a*x*x + b*x + c


def thirdorder_polynomial(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x +d


def threevariable_fit(X, a,b,c,d):
    '''
    log(R23) = ax^2 +bx +c +dlog(O32)
    '''
    x, lO32 = X
    return a * x**3 + b*x**2 + c * x + d*lO32


def LAC_two_variable(lR23_ini, lO32_ini, OH_ini, third_order=True):
    lR23 = []
    lO32 = []
    OH = []
    for ii in range(len(lR23_ini)):
        lR23 = np.concatenate([lR23, lR23_ini[ii]])
        lO32 = np.concatenate([lO32, lO32_ini[ii]])
        OH = np.concatenate([OH, OH_ini[ii]])

    OH_range = np.linspace(np.min(OH), np.max(OH), len(lR23))
    #para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)),
    # (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))

    if third_order:
        p0 = [-0.32293, 7.2954, -54.8284, 138.0430]
        o1, o2 = curve_fit(thirdorder_polynomial, OH, lR23, p0=p0)
    else:
        p0 = [7.2954, -54.8284, 138.0430]
        o1, o2 = curve_fit(secondorder_polynomial, OH, lR23, p0=p0)

    return o1, o2, lR23, lO32, OH, OH_range


def LAC_three_variable(lR23_ini, lO32_ini, OH_ini):
    '''
    log(R23) = ax^2 + bx + c + dlog(O32)
    '''
    lR23 = []
    lO32 = []
    OH = []
    for ii in range(len(lR23_ini)):
        lR23 = np.concatenate([lR23, lR23_ini[ii]])
        lO32 = np.concatenate([lO32, lO32_ini[ii]])
        OH = np.concatenate([OH, OH_ini[ii]])

    OH_range = np.linspace(np.min(OH), np.max(OH), len(lR23))
    #para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)),
    # (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))

    p0 = [7.2954, -54.8284, 138.0430, 0]
    o1, o2 = curve_fit(threevariable_fit, (OH, lO32), lR23, p0=p0)

    return o1, o2, lR23, lO32, OH, OH_range


def plot_differences_curvefit(lR23, lO32, OH, pdf_file,
                              new_coefficients=[], data_err=[], OH_err=[],
                              OH_range=[], data_range=[], marker=[], label=[],
                              IDs=[], log=None):
    """
    Plot differences between LACR23 vs observed R23 as
    a function of metallicity
    Plot differences between LACO32 vs observed O32 as
    a function of metallicity
    Used by main() function
    """
    print('new coeff: ', new_coefficients, type(new_coefficients))
    if log is None:
        log = log_stdout()

    log.info("starting ...")

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    # Plotting
    ctype = ['blue', 'green', 'red', 'magenta']

    if len(marker) == 0:
        marker = ['o'] * n_sample

    diff0 = []
    for nn in range(n_sample):
        if len(new_coefficients) != 0:
            LAC = threevariable_fit((OH[nn], lO32[nn]), *new_coefficients)
            log.info(f"curve fit LAC_R23: {LAC}")

        if nn == 0:
            if IDs:
                for jj in range(len(LAC)):
                    id_diff = lR23[nn][jj] - LAC[jj]
                    ax.annotate(IDs[nn][jj], (OH[nn][jj], id_diff),
                                fontsize='6')
        '''if n_sample == 4:
            if IDs:
                print('I got here!')
                label_IDs = [IDs[0], IDs[1]]
                for aa in range(len(label_IDs)):
                    for jj in range(len(IDs[aa])):
                        id_diff = lR23[aa][jj]-LAC[jj]
                        ax.annotate(IDs[aa][jj],
                                    (OH[aa][jj], id_diff),
                                    fontsize='6')

        else:
            if IDs:
                label_IDs = IDs[0]
                for jj in range(len(label_IDs)):
                    id_diff = lR23[0][jj] - LAC[jj]
                    ax.annotate(label_IDs[jj], (lR23[0][jj], id_diff),
                                fontsize='6')'''

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
        for ii in range(len(lR23)):
            y_ii_min = np.min(lR23[nn])
            y_ii_max = np.max(lR23[nn])
            idx = np.where((lR23[nn] >= y_ii_min) &
                           (lR23[nn] <= y_ii_max))[0]

            if len(idx) > 0:
                i_diff = lR23[nn][idx] - LAC[idx]
                ax.scatter(OH[nn], i_diff, color=ctype[nn], marker=marker[nn],
                           edgecolor='none', alpha=0.5)
                diff0 += list(lR23[nn][idx] - LAC[idx])

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


    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.set_ylabel(r'$\Delta_{R_{23}} \equiv \log(R_{23}) '
                  r'- \log(R_{23})_{\rm LAC}$')

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


def run_experiment_LAC(fitspath, fitspath_ini, secondorder=True,
                       threevariable=True, raw=False, apply_dust=False,
                       revised=True, include_rlimit=False):
    """

    Parameters
    ----------
    fitspath
    fitspath_ini
    secondorder -> used when just fitting using lR23 and OH; determines if we
                    have a second order or third order fit
    threevariable -> means that we are fitting using lR23, lO32, and OH for fit
    raw
    apply_dust
    revised
    include_rlimit

    Returns
    -------

    """
    suffix = ''
    if not revised:
        suffix += '.valid1'

    if not raw:
        suffix += '.MC'

    if apply_dust:
        suffix += '.dustcorr'

    # Validation Table Call
    if revised:
        bin_valid_file = join(fitspath, filename_dict['bin_valid_rev'])
    else:
        bin_valid_file = join(fitspath, filename_dict['bin_valid'])
    bin_valid_tab = asc.read(bin_valid_file)

    bin_derived_prop_file = join(fitspath,
                                 f"bin_derived_properties{suffix}.tbl")
    bin_derived_prop_tab = asc.read(bin_derived_prop_file)

    detect = bin_valid_tab['Detection']
    det_4363 = np.where(detect == 1)[0]
    rlimit = np.where(detect == 0.5)[0]

    # Tables of individual detections from DEEP2 and MACT samples
    derived_deep_file = join(fitspath_ini, 'DEEP2_Commons/Catalogs/',
                             'DEEP2_R23_O32_derived.tbl')
    derived_deep_tab = asc.read(derived_deep_file)

    derived_mact_file = join(fitspath_ini, 'MACT_Commons/Catalogs/',
                             'MACT_R23_O32_derived.tbl')
    derived_mact_tab = asc.read(derived_mact_file)

    # DEEP2 Derived
    # IDs are not as usefully so passing in a blank array
    DEEP2_id = []  # derived_deep_tab['ID'].data
    DEEP2_lR23 = np.log10(derived_deep_tab['R23'].data)
    DEEP2_lO32 = np.log10(derived_deep_tab['O32'].data)
    DEEP2_OH = derived_deep_tab['OH'].data

    # MACT Derived
    # IDs are not as usefully so passing in a blank array
    MACT_ID = []  # derived_mact_tab['ID'].data
    MACT_lR23 = np.log10(derived_mact_tab['R23'].data)
    MACT_lO32 = np.log10(derived_mact_tab['O32'].data)
    MACT_OH = derived_mact_tab['OH'].data

    ID = bin_derived_prop_tab['bin_ID']
    lR23_all = bin_derived_prop_tab['logR23_composite']
    lO32_all = bin_derived_prop_tab['logO32_composite']
    com_O_log = bin_derived_prop_tab['12+log(O/H)']

    det_ID = ID[det_4363]
    det_lR23 = lR23_all[det_4363]
    det_lO32 = lO32_all[det_4363]
    det_OH = com_O_log[det_4363]

    rlimit_ID = ID[rlimit]
    rlimit_lR23 = lR23_all[rlimit]
    rlimit_lO32 = lO32_all[rlimit]
    rlimit_OH = com_O_log[rlimit]

    label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']
    marker = ['D', r'$\uparrow$', '3', '4']
    color = ['b', 'g', 'r', 'm']
    #marker_a = ['D'] * len(det_lR23)
    #marker_b = [r'$\uparrow$'] * len(rlimit_lR23)
    #marker_c = ['3'] * len(DEEP2_lR23)
    #marker_d = ['4'] * len(MACT_lR23)

    #color_a = ['b'] * len(det_lR23)
    #color_b = ['g']*len(rlimit_lR23)
    #color_c = ['red'] * len(DEEP2_lR23)
    #color_d = ['magenta']*len(MACT_lR23)

    if include_rlimit:
        # color = [color_a, color_b, color_c, color_d]
        # marker = [marker_a, marker_b, marker_c, marker_d]
        lR23_arrs = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, rlimit_ID, DEEP2_id, MACT_ID]
        pdf_file = join(fitspath, f"LAC_curvefit_include_rlimit{suffix}.pdf")

    else:
        # color = [color_a, color_c, color_d]
        # marker = [marker_a, marker_c, marker_d]
        lR23_arrs = [det_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, DEEP2_id, MACT_ID]
        pdf_file = join(fitspath, f"LAC_curvefit{suffix}.pdf")

    if threevariable:
        o1, o2, lR23, lO32, OH, OH_range = LAC_three_variable(lR23_arrs, lO32_arrs, OH_arrs)
        lO32_median = np.median(lO32)
        lo32_values = [.25 * lO32_median, .75 * lO32_median, lO32_median,
                       1.25 * lO32_median, 1.75 * lO32_median]
        lo32_labels = ['lO32 = .25 * lO32_median', 'lO32 = .75 * lO32_median',
                       'lO32 = lO32_median', 'lO32 = 1.25 * lO32_median',
                       'lO32 = 1.75 * lO32_median']
        fitted_poly = np.zeros((len(lo32_values), len(lR23)))
        for aa in range(len(lo32_values)):
            fitted_poly[aa] = threevariable_fit((OH_range, lo32_values[aa]), *o1)
    else:

        if secondorder:
            o1, o2, lR23, lO32, OH, OH_range = LAC_two_variable(lR23_arrs,
                                                                lO32_arrs,
                                                                OH_arrs, third_order=False)
            fitted_poly = secondorder_polynomial(OH_range, *o1)
        else:
            o1, o2, lR23, lO32, OH, OH_range = LAC_two_variable(lR23_arrs,
                                                                lO32_arrs,
                                                                OH_arrs,
                                                                third_order=True)
            fitted_poly = thirdorder_polynomial(OH_range, *o1)

    # This is for getting the original bian plot line
    bian_R23 = local_analog_calibration.bian18_R23_OH(OH_range, bian_coeff)

    fig, ax = plt.subplots()
    # First we plot the new fitted curves and annotate with equation
    if threevariable:
        for aa in range(len(lo32_values)):
            ax.plot(fitted_poly[aa], OH_range, label=lo32_labels[aa])
        txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^2 + {o1[1]:.3f}*x"
        txt0 += rf"+ {o1[2]:.3f} + {o1[3]:.3f}*log(O32)$"
    else:
        ax.plot(fitted_poly, OH_range, label='Curve Fit')
        txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^3 + {o1[1]:.3f}*x^2"
        txt0 += rf"+ {o1[2]:.3f}*x + {o1[3]:.3f}$"
    txt0 += f"\n x = 12+log(O/H)"
    ax.annotate(txt0, [0.05, 0.92], xycoords='axes fraction', va='top',
                fontsize=6)

    # Next we plot the bian calibration
    ax.plot(bian_R23, OH_range, 'k--', label='Bian+(2018)')
    avail_idx = np.where((OH_range >= 7.80) & (OH_range <= 8.4))[0]
    ax.plot(bian_R23[avail_idx], OH_range[avail_idx], 'k-',
            label='Bian Data Range')

    # Then we plot the data
    for nn in range(len(lR23_arrs)):
        ax.scatter(lR23_arrs[nn], OH_arrs[nn],
                   marker=marker[nn], color=color[nn])

    # Next we plot the IDS
    if include_rlimit:
        # Detections and Robust limits should always be the first two arrays
        # in the list of arrays
        label_IDs = [ID_arrs[0], ID_arrs[1]]
        for aa in range(len(label_IDs)):
            for jj in range(len(ID_arrs[aa])):
                ax.annotate(ID_arrs[aa][jj],
                            (lR23_arrs[aa][jj], OH_arrs[aa][jj]), fontsize='6')
    else:
        for jj in range(len(ID_arrs[0])):
            ax.annotate(ID_arrs[0][jj], (lR23_arrs[0][jj], OH_arrs[0][jj]),
                        fontsize='6')

    ax.set(xlim=(0.0, 1.2))
    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})_{T_e}$')
    ax.legend(loc='lower left', framealpha=0.5, fontsize=6)
    fig.savefig(pdf_file)

    # This section plots makes the plot difference files
    if threevariable:
        if include_rlimit:
            if secondorder:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_threeparam_secondorder"
                                         f"_RL_diff{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_threeparam"
                                         f"_RL_diff{suffix}.pdf")
        else:
            if secondorder:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_threeparam_secondorder"
                                         f"_diff{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_threeparam_"
                                         f"_diff{suffix}.pdf")
        # print(o1[0])
        plot_differences_curvefit(lR23_arrs, lO32_arrs, OH_arrs,
                                  R23_diff_pdf_file,
                                  new_coefficients=o1, data_err=[],
                                  OH_err=[], OH_range=[np.min(OH_range),
                                                       np.max(OH_range)],
                                  data_range=[-0.5, 0.5], marker=marker,
                                  label=label, IDs=ID_arrs, log=None)
    else:
        # This is  when only doing a fit with lR23 and OH
        if include_rlimit:
            if secondorder:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_twoparam_secondorder"
                                         f"_RL_diff{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_twoparam"
                                         f"_RL_diff{suffix}.pdf")
        else:
            if secondorder:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_twoparam_secondorder"
                                         f"_diff{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath,
                                         f"LAC_curvefit_twoparam_"
                                         f"_diff{suffix}.pdf")
        local_analog_calibration.plot_differences(lR23_arrs, OH_arrs,
                                                  R23_diff_pdf_file,
                                                  data_input='R23',
                                                  new_coefficients=o1,
                                                  data_err=[],
                                                  OH_err=[],
                                                  OH_range=[np.min(OH_range),
                                                            np.max(OH_range)],
                                                  data_range=[-0.5, 0.5],
                                                  marker=marker, label=label,
                                                  IDs=[], log=None)
        local_analog_calibration.plot_differences(lO32_arrs, OH_arrs,
                                                  R23_diff_pdf_file,
                                                  data_input='O32',
                                                  new_coefficients=o1,
                                                  data_err=[],
                                                  OH_err=[],
                                                  OH_range=[np.min(OH_range),
                                                            np.max(OH_range)],
                                                  data_range=[-0.5, 0.5],
                                                  marker=marker, label=label,
                                                  IDs=[], log=None)

