import numpy as np
from astropy.io import ascii as asc
from scipy.optimize import curve_fit
from os.path import join
import matplotlib.pyplot as plt
from Zcalbase_gal.analysis import local_analog_calibration
from Metallicity_Stack_Commons.column_names import filename_dict
from Zcalbase_gal.analysis.deep2_r23_o32 import bian_coeff, \
    secondorder_polynomial, thirdorder_polynomial, threevariable_fit
from .plotting import plot_difference_curvefit
from . import ctype
from .log_commons import log_stdout, LogClass


def fitting_function(lR23_array, lO32_array, OH_array, secondorder=True,
                     threevariable=True):
    """
    :param lR23_array: array. an array of arrays of log(R23) values; each array
                        from a different data set
    :param lO32_array: array. an array of arrays of log(O32) values; each array
                        from a different data set
    :param OH_array: array. an array of arrays of OH values; each array
                        from a different data set
    Returns
    o1, o2 -> contains all parameters from teh curve fitting
    lR23, lO32, OH -> array of all specific values concatenated and no longer
                        distinguishable by data sete
    OH_range -> array of OH values that span the range of OH
    """
    lR23 = []
    lO32 = []
    OH = []
    for ii in range(len(lR23_array)):
        lR23 = np.concatenate([lR23, lR23_array[ii]])
        lO32 = np.concatenate([lO32, lO32_array[ii]])
        OH = np.concatenate([OH, OH_array[ii]])
    OH_range = np.linspace(np.min(OH), np.max(OH), len(lR23))

    # Currently not using a threevariable third order fit which can be added
    if threevariable:
        p0 = [7.2954, -54.8284, 138.0430, 0]
        o1, o2 = curve_fit(threevariable_fit, (OH, lO32), lR23, p0=p0)
    else:
        if secondorder:
            p0 = [7.2954, -54.8284, 138.0430]
            o1, o2 = curve_fit(secondorder_polynomial, OH, lR23, p0=p0)
        else:
            p0 = [-0.32293, 7.2954, -54.8284, 138.0430]
            o1, o2 = curve_fit(thirdorder_polynomial, OH, lR23, p0=p0)

    return o1, o2, lR23, lO32, OH, OH_range


def run_experiment_Zcal(fitspath, fitspath_curvefit, fitspath_ini, n_bins=4,
                        secondorder=True, threevariable=True, raw=False,
                        apply_dust=False, revised=True, include_rlimit=False):
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

    threevariable = True, secondorder = True
    lR23 = ax^2+ bx+ c+ dlO32
    threevariable = True, secondorder = False
    lR23 = gx^3 + ax^2+ bx+ c+ dlO32
    threevariable = False, secondorder = True
    lR23 = ax^2+ bx+ c
    threevariable = False, secondorder = False
    lR23 = gx^3 + ax^2+ bx+ c
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

    # Names of files
    if threevariable:
        if include_rlimit:
            if secondorder:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"threeparam_secondorder"
                                         f"_includerobust_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"threeparam_secondorder_includerobust"
                                f"{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"threeparam_thirdorder"
                                         f"_includerobust_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"threeparam_thirdorder_includerobust{suffix}.pdf")
        else:
            if secondorder:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"threeparam_secondorder"
                                         f"_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"threeparam_secondorder{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"threeparam_thirdorder"
                                         f"_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"threeparam_thirdorder{suffix}.pdf")
    else:
        # This is  when only doing a fit with lR23 and OH
        if include_rlimit:
            if secondorder:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"twoparam_secondorder"
                                         f"_includerobust_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"twoparam_secondorder_RL_"
                                f"{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"twoparam_thirdorder"
                                         f"_includerobust_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"twoparam_thirdorder_includerobust{suffix}.pdf")
        else:
            if secondorder:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"twoparam_secondorder"
                                         f"_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"twoparam_secondorder{suffix}.pdf")
            else:
                R23_diff_pdf_file = join(fitspath_curvefit,
                                         f"twoparam_thirdorder"
                                         f"_diff{suffix}.pdf")
                pdf_file = join(fitspath_curvefit,
                                f"twoparam_thirdorder{suffix}.pdf")

    if include_rlimit:
        lR23_arrs = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, rlimit_ID, DEEP2_id, MACT_ID]
    else:
        lR23_arrs = [det_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, DEEP2_id, MACT_ID]

    if threevariable:
        o1, o2, lR23, lO32, OH, OH_range = \
            fitting_function(lR23_arrs, lO32_arrs, OH_arrs, secondorder=True,
                             threevariable=True)
    else:

        if secondorder:
            o1, o2, lR23, lO32, OH, OH_range = \
                fitting_function(lR23_arrs, lO32_arrs, OH_arrs,
                                 secondorder=True, threevariable=False)
            fitted_poly = secondorder_polynomial(OH_range, *o1)
        else:
            o1, o2, lR23, lO32, OH, OH_range = \
                fitting_function(lR23_arrs, lO32_arrs, OH_arrs,
                                 secondorder=False, threevariable=False)
            fitted_poly = thirdorder_polynomial(OH_range, *o1)

    # Creating the inning by lO32 values for plots
    bin_start = np.zeros(n_bins)
    bin_end = np.zeros(n_bins)

    sort0 = np.argsort(lO32)
    y_sort0 = lO32[sort0]
    r_bin_pts = np.int(np.round(len(lO32) / float(n_bins)))

    lo32_bin_avg = np.zeros(n_bins)
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
        lo32_bin_avg[ii] = np.average(y_sort0[lo32_in_bin])
    print(lo32_bin_avg)
    fitted_poly = np.zeros((len(lo32_bin_avg), len(lR23)))
    print('*o1', *o1)
    if threevariable:
        # ax^2 + bx +c + dlog(O32)
        for aa in range(len(lo32_bin_avg)):
            fitted_poly[aa] = threevariable_fit((OH_range, lo32_bin_avg[aa]),
                                                *o1)
    else:
        if secondorder:
            # ax^2+bx+c
            fitted_poly = secondorder_polynomial(OH_range, *o1)
        else:
            # ax^3 + bx^2 +cx+d
            fitted_poly = thirdorder_polynomial(OH_range, *o1)

    x_range = [-0.5, 0.5]
    # This is for getting the original bian plot line
    bian_R23 = local_analog_calibration.bian18_R23_OH(OH_range, bian_coeff)

    # Starting the plotting
    fig, ax = plt.subplots()

    # Plotting the binned data and its legend
    for nn in range(len(lR23_arrs)):
        for ii in range(n_bins):
            y_ii_min = bin_start[ii]  # bin_y_min + ii * dy
            y_ii_max = bin_end[ii]  # y_min + (ii+1) * dy
            # print('bin start: ', bin_start[ii], 'bin end: ', bin_end[ii])
            idx = np.where((lO32_arrs[nn] >= y_ii_min) &
                           (lO32_arrs[nn] <= y_ii_max))[0]
            ii_label = ''
            if nn == len(lR23_arrs)-1:
                idx_all = np.where((lO32 >= y_ii_min) & (lO32 <= y_ii_max))[0]
                ii_label = fr" {y_ii_min:.2f} < $\log(O_{{32}})$ " + \
                           f"< {y_ii_max:.2f}, N = {len(idx_all):d}"
            if len(idx) > 0:
                ax.scatter(lR23_arrs[nn][idx], OH_arrs[nn][idx],
                           color=ctype[ii], marker=marker[nn], label=ii_label)

    # Now we plot the new fitted curves and annotate with equation
    if threevariable:
        # ax^2 + bx +c + dlog(O32)
        for aa in range(len(lo32_bin_avg)):
            ax.plot(fitted_poly[aa], OH_range, color=ctype[aa],
                    label=f"Average lO32: {lo32_bin_avg[aa]:.3f}")
        txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^2 + {o1[1]:.3f}*x"
        txt0 += rf"+ {o1[2]:.3f} + {o1[3]:.3f}*log(O32)$"
    else:
        if secondorder:
            # ax^2+bx+c
            ax.plot(fitted_poly, OH_range, color=ctype[-1], label='Curve Fit')
            txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^2 + {o1[1]:.3f}*x"
            txt0 += rf"+ {o1[2]:.3f}$"
        else:
            # ax^3 + bx^2 +cx+d
            ax.plot(fitted_poly, OH_range, color=ctype[-1], label='Curve Fit')
            txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^3 + {o1[1]:.3f}*x^2"
            txt0 += rf"+ {o1[2]:.3f}*x + {o1[3]:.3f}$"

    txt0 += f"\n x = 12+log(O/H)"
    ax.annotate(txt0, [0.05, 0.92], xycoords='axes fraction', va='top',
                fontsize=6)
    ax.legend(loc='lower left', framealpha=0.5, fontsize=6)

    # Next we plot the bian calibration
    ax.plot(bian_R23, OH_range, 'k--', label='Bian+(2018)')
    avail_idx = np.where((OH_range >= 7.80) & (OH_range <= 8.4))[0]
    ax.plot(bian_R23[avail_idx], OH_range[avail_idx], 'k-',
            label='Bian Data Range')

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
    # print('PDF file name: ', pdf_file)
    fig.savefig(pdf_file)
    fitting_model = 'Zcal'
    # This section plots makes the plot difference files
    if threevariable:
        plot_difference_curvefit.\
            plot_difference_threevariable(lR23_arrs, lO32_arrs, OH_arrs, lO32,
                                          R23_diff_pdf_file, fitting_model,
                                          bin_start, bin_end,
                                          new_coefficients=o1, n_bins=4,
                                          data_err=[], OH_err=[],
                                          OH_range=[np.min(OH_range),
                                                    np.max(OH_range)],
                                          data_range=[-0.3, 0.3],
                                          marker=marker, label=label,
                                          IDs=ID_arrs, log=None)
    else:
        plot_difference_curvefit.\
            plot_difference_twovariable(lR23_arrs, lO32_arrs, OH_arrs,
                                        bin_start, bin_end, R23_diff_pdf_file,
                                        new_coefficients=o1, n_bins=4,
                                        data_err=[], OH_err=[],
                                        OH_range=[np.min(OH_range),
                                                  np.max(OH_range)],
                                        data_range=[-0.3, 0.3],
                                        marker=marker, label=label, IDs=[],
                                        log=None)
