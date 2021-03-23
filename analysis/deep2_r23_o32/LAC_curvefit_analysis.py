import numpy as np
from astropy.io import ascii as asc
from scipy.optimize import curve_fit
from os.path import join
import matplotlib.pyplot as plt
from Zcalbase_gal.analysis import local_analog_calibration
from Metallicity_Stack_Commons.column_names import filename_dict
from Zcalbase_gal.analysis.deep2_r23_o32 import bian_coeff


def secondorder_polynomial(x, a, b, c):
    return a*x*x + b*x + c


def thirdorder_polynomial(x, a, b, c, d):
    return a*x*x*x + b*x*x + c*x +d


def threevariable_fit(X, a, b, c, d):
    x, lO32 = X
    return a * x * x * x + b * x * x + c * x + d*lO32


def LAC_two_variable(lR23_ini, lO32_ini, OH_ini, order):
    lR23 = []
    lO32 = []
    OH = []
    for ii in range(len(lR23_ini)):
        lR23 = np.concatenate([lR23, lR23_ini[ii]])
        lO32 = np.concatenate([lO32, lO32_ini[ii]])
        OH = np.concatenate([OH, OH_ini[ii]])

    OH_range = np.linspace(np.min(OH), np.max(OH), 100)

    #para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)),
    # (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))

    if order == 'second':
        p0 = [-0.32293, 7.2954, -54.8284, 138.0430]
        o1, o2 = curve_fit(thirdorder_polynomial, OH, lR23, p0=p0)
    else:
        p0 = [7.2954, -54.8284, 138.0430]
        o1, o2 = curve_fit(thirdorder_polynomial, OH, lR23, p0=p0)

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

    OH_range = np.linspace(np.min(OH), np.max(OH), 100)

    #para_bound = ((working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0)),
    # (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0)))

    p0 = [7.2954, -54.8284, 138.0430, 0]
    o1, o2 = curve_fit(threevariable_fit, (OH, lO32), lR23, p0=p0)

    return o1, o2, lR23, lO32, OH, OH_range


def run_experiment_LAC(fitspath, fitspath_ini, secondorder=True,
                       threevariable=True, raw=False, apply_dust=False,
                       revised=False, include_rlimit=False):
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

    marker_a = ['D'] * len(det_lR23)
    marker_b = [r'$\uparrow$'] * len(rlimit_lR23)
    marker_c = ['3'] * len(DEEP2_lR23)
    marker_d = ['4'] * len(MACT_lR23)

    color_a = ['b'] *len(det_lR23)
    color_b = ['g']*len(rlimit_lR23)
    color_c = ['red'] * len(DEEP2_lR23)
    color_d = ['magenta']*len(MACT_lR23)

    if include_rlimit:
        color = np.concatenate([color_a, color_b, color_c, color_d])
        marker = np.concatenate([marker_a, marker_b, marker_c, marker_d])
        lR23_arrs = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, rlimit_ID, DEEP2_id, MACT_ID]
        pdf_file = join(fitspath, f"LAC_curvefit_include_rlimit{suffix}.pdf")
        R23_diff_pdf_file = join(fitspath,
                                 f"LAC_curvefit_R23_include_rlimit_diff{suffix}.pdf")
    else:
        color = np.concatenate([color_a, color_c, color_d])
        marker = np.concatenate([marker_a, marker_c, marker_d])
        lR23_arrs = [det_lR23, DEEP2_lR23, MACT_lR23]
        lO32_arrs = [det_lO32, DEEP2_lO32, MACT_lO32]
        OH_arrs = [det_OH, DEEP2_OH, MACT_OH]
        ID_arrs = [det_ID, DEEP2_id, MACT_ID]
        pdf_file = join(fitspath, f"LAC_curvefit{suffix}.pdf")
        R23_diff_pdf_file = join(fitspath, f"LAC_curvefit_R23_diff{suffix}.pdf")

    if threevariable:
        o1, o2, lR23, lO32, OH, OH_range = LAC_three_variable(lR23_arrs,
                                                              lO32_arrs,
                                                              OH_arrs)
    else:
        o1, o2, lR23, lO32, OH, OH_range = LAC_two_variable(lR23_arrs, lO32_arrs,
                                                            OH_arrs)
    print('o1: ', o1)
    if secondorder == 'second':
        fitted_poly = secondorder_polynomial(OH_range, *o1)
        bian_R23 = local_analog_calibration.bian18_R23_OH(OH_range, bian_coeff)
    else:
        fitted_poly = thirdorder_polynomial(OH_range, *o1)
        bian_R23 = local_analog_calibration.bian18_R23_OH(OH_range, bian_coeff)

    fig, ax = plt.subplots()
    ax.plot(fitted_poly, OH_range, label='Curve Fit')
    ax.plot(bian_R23, OH_range, 'k--', label='Bian+(2018)')
    avail_idx = np.where((OH_range >= 7.80) & (OH_range <= 8.4))[0]
    ax.plot(bian_R23[avail_idx], OH_range[avail_idx], 'k-',
            label='Bian Data Range')
    for ii in range(len(marker)):
        ax.scatter(lR23[ii], OH[ii], marker=marker[ii], color=color[ii])
    if threevariable:
        txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^2 + {o1[1]:.3f}*x"
        txt0 += rf"+ {o1[2]:.3f} + {o1[3]:.3f}*log(O32)$"
    else:
        txt0 = rf"curve_fit: $log(R23) = {o1[0]:.3f}*x^3 + {o1[1]:.3f}*x^2"
        txt0 += rf"+ {o1[2]:.3f}*x + {o1[3]:.3f}$"
    txt0 += f"\n x = 12+log(O/H)"
    ax.annotate(txt0, [0.05, 0.92], xycoords='axes fraction', va='top', fontsize=6)
    ax.set(xlim=(0.0, 1.2))
    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})_{T_e}$')
    ax.legend(loc='lower left', framealpha=0.5, fontsize=10)

    local_analog_calibration.plot_differences(lR23_arrs, OH_arrs, R23_diff_pdf_file,
                                              data_input='R23',
                                              new_coefficients=o1, data_err=[],
                                              OH_err=[], OH_range=[np.min(OH_range), np.max(OH_range)],
                                              data_range=[-0.5, 0.5],
                                              marker=marker, label=label,
                                              IDs=[], log=None)
    fig.savefig(pdf_file)
