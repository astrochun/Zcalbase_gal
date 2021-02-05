# Graphs the temperature, metallicities,
# R23 and O32 and errors for the individual and
# composite spectra by importing pre-existing tables and dictionaries

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from ..log_commons import log_stdout

from Metallicity_Stack_Commons.column_names import filename_dict

jiang_coeff = [-24.135, 6.1532, -0.37866, -0.147, -7.071]
bian_coeff = [138.0430, -54.8284, 7.2954, -0.32293]


def jiang_calibration(metal_det, lo32):
    """
    Calculating the Jiang calibration for the metallicity
    and O32 measurements enter
    Jiang Calibration:
    jiang_coeff[0] + jiang_coeff[1] * metal_det + jiang_coeff[2]
    * (metal_det * metal_det) + jiang_coeff[3]
    * (jiang_coeff[4] + metal_det) * lo32

    :param metal_det: dict. dictionary of 12+log(O/H)
    :param lO32: np.array. log(O32) value

    Returns Jiang's expected R23 value given metallicity and O32
    """

    p_jiang = np.poly1d([jiang_coeff[2], jiang_coeff[1], jiang_coeff[0]])
    shortcut = p_jiang(metal_det) + \
               jiang_coeff[3] * (jiang_coeff[4] + metal_det) * lo32
    return shortcut


def bian_calibration_r23(metal_det):
    """
    Calculating the R23 of Bian calibration for the
    metallicity measurements enter
    Bian R23 Calibration:
    bian_coeff[0] + bian_coeff[1] * metal_det + bian_coeff[2]
            *metal_det*metal_det+ bian_coeff[3]
            * metal_det * metal_det * metal_det

    :param metal_det: dict. dictionary of 12+log(O/H)

    Returns Bian's expected R23 value given metallicity
    """
    p_bian = np.poly1d([bian_coeff[::-1]])
    shortcut = p_bian(metal_det)
    return shortcut


def bian_calibration_o32(metal_det):
    """
    Calculating the O32 of Bian calibration for
    the metallicity measurements enter

    :param metal_det: dict. dictionary of 12+log(O/H)

    Returns Bian's expected O32 values given metallicity
    """
    return (-1 / 0.59) * (metal_det - 8.54)


def plotting_te_metal(fitspath, fitspath_ini, raw=False, apply_dust=False,
                      revised=False, individual=False, log=None):
    """
    Plotting metallicity, temperature, R23, and O32 for different datasets

    :param fitspath: str. path where files are called from and saved to
    :param fitspath_ini: str. path where whole project is stored
    :param raw: bool.
    :param apply_dust: bool.
    :param revised: bool. keyword (automatically set to False) to determine
                        what data is used refers to if using the
                        bin_derived_prop_revised temperature and metallicity
                        measurements which right now implement dust attenuation
    :param individual: bool.
    :param log: LogClass. Default use log_stdout()

    PDF File: fitspath + temperature_metallicity_plots{suffix}.pdf"
    No returns
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    suffix = ''
    if not revised:
        suffix += '.valid1'

    if not raw:
        suffix += '.MC'

    if apply_dust:
        suffix += '.dustcorr'

    out_pdf = join(fitspath, f"temperature_metallicity_plots{suffix}.pdf")

    temperature_table = join(fitspath, f"bin_derived_properties{suffix}.tbl")
    log.info(f"Reading: {temperature_table}")
    comp_derived = asc.read(temperature_table)

    if revised:
        verification = asc.read(join(fitspath, filename_dict['bin_valid_rev']))
    else:
        verification = asc.read(join(fitspath, filename_dict['bin_valid']))
    log.info(f"Reading: {verification}")

    # Individual Measurements
    if individual:
        indv_all_file = join(fitspath, filename_dict['indv_derived_prop'])
        log.info(f"Reading: {indv_all_file}")
        indv_all = asc.read(indv_all_file)

        icom_log = indv_all['12+log(O/H)']
        ilogR23 = indv_all['logR23']
        ilogO32 = indv_all['logO32']
        iidx = np.where((icom_log != 0.0))[0]

        iiR23_idv = ilogR23[iidx]
        iiO32_idv = ilogO32[iidx]
        iR23_idv = np.log10(iiR23_idv)
        iO32_idv = np.log10(iiO32_idv)
        log.info(f"len: {len(iR23_idv)}")

    # DEEP2 and MACT Data
    DEEP_table = join(fitspath_ini,
                      'DEEP2_Commons/Catalogs/DEEP2_R23_O32_derived.tbl')
    MACT_table = join(fitspath_ini,
                      'MACT_Commons/Catalogs/MACT_R23_O32_derived.tbl')
    log.info(f"Reading: {DEEP_table}")
    derived_DEEP = asc.read(DEEP_table)

    log.info(f"Reading: {MACT_table}")
    derived_MACT = asc.read(MACT_table)

    # DEEP2 Derived
    der_R23 = np.log10(derived_DEEP['R23'].data)
    der_O32 = np.log10(derived_DEEP['O32'].data)
    der_Te = derived_DEEP['Te'].data
    der_OH = derived_DEEP['OH'].data
    ID_der = derived_DEEP['ID'].data
    # der_OH_log = np.log10(er_OH_log)

    # MACT Derived
    der_R23_MACT = np.log10(derived_MACT['R23'].data)
    der_O32_MACT = np.log10(derived_MACT['O32'].data)
    der_Te_MACT = derived_MACT['Te'].data
    der_OH_MACT = derived_MACT['OH'].data
    ID_der_MACT = derived_MACT['ID'].data

    # Composite Measurements
    R23_composite = comp_derived['logR23_composite'].data
    O32_composite = comp_derived['logO32_composite'].data
    ID_composite = comp_derived['bin_ID'].data
    T_e_composite = comp_derived['T_e'].data
    metal_composite = comp_derived['12+log(O/H)'].data

    # Verification Table
    ver_detection = verification['Detection'].data
    ver_detect = np.where((ver_detection == 1))[0]
    ver_rlimit = np.where((ver_detection == 0.5))[0]

    pdf_pages = PdfPages(out_pdf)
    if individual:
        fig, ax = plt.subplots()
        ax.scatter(iR23_idv, iO32_idv, marker='.', s=35, color='g')
        ax.set_title(r'$R_{23}$ vs. $O_{32}$')
        ax.set_xlabel(r'log($R_{23}$)')
        ax.set_ylabel(r'log($O_{32}$)')
        fig.savefig(pdf_pages, format='pdf')
        fig.clear()

    # ########################################################################
    fig1, ax1 = plt.subplots()

    ax1.scatter(T_e_composite[ver_detect], R23_composite[ver_detect],
                marker='.', s=50, color='b')
    ax1.scatter(T_e_composite[ver_rlimit], R23_composite[ver_rlimit],
                marker='<', s=35, color='b')
    for xx in ver_detect:
        ax1.annotate(ID_composite[xx], (T_e_composite[xx], R23_composite[xx]),
                     fontsize='6')
    for xx in ver_rlimit:
        ax1.annotate(ID_composite[xx], (T_e_composite[xx], R23_composite[xx]),
                     fontsize='6')

    ax1.scatter(der_Te, der_R23, s=20, marker='*', color='k',
                edgecolors='None')
    for b in range(len(ID_der)):
        ax1.annotate(ID_der[b], (der_Te[b], der_R23[b]), fontsize='2')

    ax1.scatter(der_Te_MACT, der_R23_MACT, s=20, marker='P', color='r',
                alpha=0.5, edgecolors='None')
    for q in range(len(ID_der_MACT)):
        ax1.annotate(ID_der_MACT[q], (der_Te_MACT[q], der_R23_MACT[q]),
                     fontsize='2')

    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel(r'$R_{23}$')
    ax1.set_title(r'Temperatures vs $R_{23}$ Temperature')

    fig1.savefig(pdf_pages, format='pdf')
    fig1.clear()

    # ########################################################################
    fig2, ax2 = plt.subplots()

    ax2.scatter(T_e_composite[ver_detect], O32_composite[ver_detect],
                marker='.', s=50, color='b')
    ax2.scatter(T_e_composite[ver_rlimit], O32_composite[ver_rlimit],
                marker='<', s=35, color='b')
    for c in ver_detect:
        ax2.annotate(ID_composite[c], (T_e_composite[c], O32_composite[c]),
                     fontsize='6')
    for c in ver_rlimit:
        ax2.annotate(ID_composite[c], (T_e_composite[c], O32_composite[c]),
                     fontsize='6')

    ax2.scatter(der_Te, der_O32, s=20, marker='*', color='k',
                edgecolors='None')
    for f in range(len(ID_der)):
        ax2.annotate(ID_der[f], (der_Te[f], der_O32[f]), fontsize='2')

    ax2.scatter(der_Te_MACT, der_O32_MACT, s=20, marker='P', color='r',
                alpha=0.5, edgecolors='None')
    for s in range(len(ID_der_MACT)):
        ax2.annotate(ID_der_MACT[s], (der_Te_MACT[s], der_O32_MACT[s]),
                     fontsize='2')

    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel(r'$O_{32}$')
    ax2.set_title(r'Temperatures vs $O_{32}$ Temperature')

    fig2.savefig(pdf_pages, format='pdf')
    fig2.clear()

    # ########################################################################
    fig3, ax3 = plt.subplots()
    ax3.scatter(R23_composite[ver_detect], metal_composite[ver_detect],
                marker='.', s=50, color='b')
    ax3.scatter(R23_composite[ver_rlimit], metal_composite[ver_rlimit],
                marker='^', s=35, color='b')
    for zz in ver_detect:
        ax3.annotate(ID_composite[zz], (R23_composite[zz],
                                        metal_composite[zz]), fontsize='6')
    for zz in ver_rlimit:
        ax3.annotate(ID_composite[zz], (R23_composite[zz],
                                        metal_composite[zz]), fontsize='6')

    ax3.scatter(der_R23, der_OH, s=20, marker='*', color='k',
                edgecolors='None')
    for gg in range(len(ID_der)):
        ax3.annotate(ID_der[gg], (der_R23[gg], der_OH[gg]), fontsize='2')
    ax3.scatter(der_R23_MACT, der_OH_MACT, s=20, marker='P', color='r',
                alpha=0.5, edgecolors='None')
    for g in range(len(ID_der_MACT)):
        ax3.annotate(ID_der_MACT[g], (der_R23_MACT[g], der_OH_MACT[g]),
                     fontsize='2')
    ax3.set_xlim(0.5, 1.1)
    ax3.set_ylim(6.75, 9.25)
    ax3.set_xlabel(r'$R_{23}$')
    ax3.set_ylabel('12+log(O/H)')
    ax3.set_title(r'$R_{23}$ vs. Composite Metallicity')

    fig3.savefig(pdf_pages, format='pdf')
    fig3.clear()

    # ########################################################################
    fig4, ax4 = plt.subplots()
    ax4.scatter(O32_composite[ver_detect], metal_composite[ver_detect],
                marker='.', s=50, color='b')
    ax4.scatter(O32_composite[ver_rlimit], metal_composite[ver_rlimit],
                marker='^', s=35, color='b')

    for ww in ver_detect:
        ax4.annotate(ID_composite[ww], (O32_composite[ww],
                                        metal_composite[ww]),
                     fontsize='6')
    for ww in ver_rlimit:
        ax4.annotate(ID_composite[ww], (O32_composite[ww],
                                        metal_composite[ww]),
                     fontsize='6')

    ax4.scatter(der_O32, der_OH, s=20, marker='*', color='k',
                edgecolors='None')
    for hh in range(len(ID_der)):
        ax4.annotate(ID_der[hh], (der_O32[hh], der_OH[hh]), fontsize='2')

    ax4.scatter(der_O32_MACT, der_OH_MACT, s=20, marker='P', color='r',
                alpha=0.5, edgecolors='None')
    for h in range(len(ID_der_MACT)):
        ax4.annotate(ID_der_MACT[h], (der_O32_MACT[h], der_OH_MACT[h]),
                     fontsize='2')

    ax4.set_xlabel(r'$O_{32}$')
    ax4.set_ylabel('12+log(O/H)')
    ax4.set_title(r'$O_{32}$ vs. Composite Metallicity')
    fig4.savefig(pdf_pages, format='pdf')
    fig4.clear()

    # ########################################################################
    log.info(f"Writing: {out_pdf}")
    pdf_pages.close()

    log.info("finished.")


def jiang_comparison(fitspath, fitspath_ini, log=None):
    """
    To calculate the difference between the predicted Jiang calibration
    and the measured values  for composite spectra
    Jiang Calibration fron Jiang, T., Malhotra, S., Rhoads,
    J. E., et al. 2019, ApJ, 872, 145
    Jiang Calibration
    log(R23) = a +bx+cx^2 - d(e+x)y
    x = 12+log(O/H)
    y = log(O32)

    :param fitspath: str. path to where files come and are saved to
    :param fitspath_ini: str. path where whole project is stored
    :param log: LogClass. Default use log_stdout()

    Still trying to decide if this is necessary: 
    revised -> keyword (automatically False) to determine what data is used

    PDF File: fitspath + 'comparsion_Jiang_Zcal.pdf'
    No returns
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    valid_file = join(fitspath, 'bin_validation.revised.tbl')
    log.info(f"Reading: {valid_file}")
    validation = asc.read(valid_file)

    temp_file = join(fitspath, 'bin_derived_properties.tbl')
    log.info(f"Reading: {temp_file}")
    temp_tab = asc.read(temp_file)

    out_pdf = join(fitspath, 'comparsion_Jiang_Zcal.pdf')
    pdf_pages = PdfPages(out_pdf)

    bin_id = temp_tab['bin_ID'].data
    lr23_comp = temp_tab['logR23_Composite'].data
    lo32_comp = temp_tab['logO32_Composite'].data
    zmetal = temp_tab['12+log(O/H)'].data

    valid = validation['Detection'].data
    detect = np.where(valid == 1.0)[0]
    rlimit = np.where(valid == 0.5)[0]

    valid_id = bin_id[detect]

    lR23 = lr23_comp[detect]
    lO32 = lo32_comp[detect]

    rlr23 = lr23_comp[rlimit]
    rlo32 = lo32_comp[rlimit]

    metal_det = zmetal[detect]
    metal_rl = zmetal[rlimit]

    derived_deep2_file = join(fitspath_ini, 'DEEP2_R23_O32_derived.tbl')
    log.info(f"Reading: {derived_deep2_file}")
    derived_deep2 = asc.read(derived_deep2_file)

    derived_mact_file = join(fitspath_ini, 'MACT_R23_O32_derived.tbl')
    log.info(f"Reading: {derived_mact_file}")
    derived_mact = asc.read(derived_mact_file)

    # DEEP2 Derived
    er_r23 = derived_deep2['R23'].data
    er_o32 = derived_deep2['O32'].data
    der_r23 = np.log10(er_r23)
    der_o32 = np.log10(er_o32)
    der_oh = derived_deep2['OH'].data

    # MACT Derived
    er_r23_mact = derived_mact['R23'].data
    er_o32_mact = derived_mact['O32'].data
    der_r23_mact = np.log10(er_r23_mact)
    der_o32_mact = np.log10(er_o32_mact)
    der_oh_mact = derived_mact['OH'].data

    jR23_det = jiang_calibration(metal_det, lO32)
    A_comparison = (jR23_det - lR23)
    count = len(metal_det) + len(der_r23) + len(der_r23_mact)
    log.info(f"count: {count}")

    jR23_DEEP = jiang_calibration(der_oh, der_o32)
    jR23_MACT = jiang_calibration(der_oh_mact, der_o32_mact)

    B_comparison = jR23_DEEP - der_r23
    C_comparison = jR23_MACT - der_r23_mact

    arr_sum = np.concatenate((A_comparison, B_comparison, C_comparison),
                             axis=None)
    log.info(f"arr_sum: {arr_sum}")
    med0 = np.median(arr_sum)
    avg0 = np.average(arr_sum)
    sig0 = np.std(arr_sum)
    log.info(f'med: {med0} avg: {avg0} sig: {sig0}')

    fig, ax = plt.subplots()
    ax.scatter(lR23, jR23_det, marker='D', color='b', alpha=0.75,
               label='Composite Detections')
    for aa in range(len(valid_id)):
        ax.annotate(valid_id[aa], (lR23[aa], jR23_det[aa]), fontsize='6')
    ax.scatter(der_r23, jR23_DEEP, marker='3', color='r',
               label='DEEP2 Individual Spectra')
    ax.scatter(der_r23_mact, jR23_MACT, marker='4', color='m',
               label='MACT Individual Spectra')
    ax.set_xlabel(r'Observed $log(R_{23})$')
    ax.set_ylabel(r'Estimated $log(R_{23})$')
    plt.plot(lR23, lR23, 'k', label='One to one line')
    plt.legend()

    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.annotate(an_txt, [0.2, 0.015], xycoords='axes fraction', va='bottom',
                ha='right', fontsize=10)

    fig.savefig(pdf_pages, format='pdf')
    log.info(f"Writing: {out_pdf}")
    pdf_pages.close()

    log.info("finished.")


def bian_comparison(fitspath, fitspath_ini, log=None):
    """
    To calculate the difference between the predicted Bian calibration
    and the measured values
    for composite spectra
    Bian Calibration from Bian, F., Kewley, L. J., & Dopita,
    M. A. 2018, ApJ, 859, 175

    :param fitspath: str. path to where files come and are saved to
    :param fitspath_ini: str. path where whole project is stored
    :param log: LogClass. Default use log_stdout()

    PDF File: fitspath + 'comparsion_Bian_Zcal.pdf'
    No Return
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    # Do we need an option for revised?
    valid_file = join(fitspath, 'bin_validation.revised.tbl')
    log.info(f"Reading: {valid_file}")
    validation = asc.read(valid_file)

    temp_file = join(fitspath, 'bin_derived_properties.tbl')
    log.info(f"Reading; {temp_file}")
    temp_tab = asc.read(temp_file)

    out_pdf = join(fitspath, 'comparsion_Bian_Zcal.pdf')
    pdf_pages = PdfPages(out_pdf)

    bin_ID = temp_tab['bin_ID']
    lR23_comp = temp_tab['logR23_Composite']
    lO32_comp = temp_tab['logO32_Composite']
    zmetal = temp_tab['12+log(O/H)']

    valid = validation['Detection']
    detect = np.where((valid == 1.0))[0]

    valid_ID = bin_ID[detect]

    lR23 = lR23_comp[detect]
    lO32 = lO32_comp[detect]

    metal_det = zmetal[detect]

    derived = asc.read(join(fitspath_ini, 'DEEP2_R23_O32_derived.tbl'))
    derived_MACT = asc.read(join(fitspath_ini, 'MACT_R23_O32_derived.tbl'))

    # DEEP2 Derived
    deep2_r23 = np.log10(derived['R23'].data)
    deep2_o32 = np.log10(derived['O32'].data)
    deep2_oh = derived['OH'].data

    # MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_MACT['OH'].data

    jR23_det = bian_calibration_r23(metal_det)
    jO32_det = bian_calibration_o32(metal_det)

    A_comparison = jR23_det - lR23
    AO_comparison = jO32_det - lO32

    bR23_DEEP = bian_calibration_r23(deep2_oh)
    bO32_DEEP = bian_calibration_o32(deep2_oh)
    bR23_MACT = bian_calibration_r23(der_OH_MACT)
    bO32_MACT = bian_calibration_o32(der_OH_MACT)

    B_comparison = bR23_DEEP - deep2_r23
    BO_comparison = bO32_DEEP - deep2_o32
    C_comparison = bR23_MACT - der_R23_MACT
    CO_comparison = bO32_MACT - der_O32_MACT

    arr_sum = np.concatenate((A_comparison, B_comparison, C_comparison),
                             axis=None)
    med0 = np.median(arr_sum)
    avg0 = np.average(arr_sum)
    sig0 = np.std(arr_sum)

    log.info(f"DEEP2 x: {deep2_r23}")
    log.info(f"DEEP2 y: {bR23_DEEP}")
    log.info(f"MACT x: {der_R23_MACT}")
    log.info(f"MACT y: {bR23_MACT}")
    log.info(f"concatenate array: {arr_sum}")
    log.info(f"med: {med0}  avg: {avg0}  sig: {sig0}")

    npz_outfile = join(fitspath, 'bian_comparison_xandy_values.npz')
    log.info(f"Writing: {npz_outfile}")
    np.savez(npz_outfile, DEEPx=deep2_r23, DEEPy=bR23_DEEP,
             MACTx=der_R23_MACT, MACTy=bR23_MACT)

    fig, ax = plt.subplots()
    ax.scatter(lR23, jR23_det, marker='D', color='b', alpha=0.75,
               label='Detections')
    ax.scatter(deep2_r23, bR23_DEEP, marker='3', color='r',
               label='DEEP2 Individual Spectra')
    ax.scatter(der_R23_MACT, bR23_MACT, marker='4', color='m',
               label='MACT Individual Spectra')
    for aa in range(len(valid_ID)):
        ax.annotate(valid_ID[aa], (lR23[aa], jR23_det[aa]), fontsize='6')
    ax.set_xlabel('Observed ' + r'$log(R_{23})$')
    ax.set_ylabel('Estimated ' + r'$log(R_{23})$')
    plt.plot(lR23, lR23, 'k', label='One to one line')
    plt.legend()

    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.annotate(an_txt, [0.2, 0.85], xycoords='axes fraction', va='bottom',
                ha='right', fontsize=10)

    fig.savefig(pdf_pages, format='pdf')
    fig.clear()

    arr_sum1 = np.concatenate((AO_comparison, BO_comparison, CO_comparison),
                              axis=None)
    medO0 = np.median(arr_sum1)
    avgO0 = np.average(arr_sum1)
    sigO0 = np.std(arr_sum1)

    log.info(f"Writing: {arr_sum1}")
    log.info(f"med: {medO0}  avg: {avgO0}  sig: {sigO0}")

    fig, ax = plt.subplots()
    ax.scatter(lO32, jO32_det, marker='D', color='b', label='Detections')
    ax.scatter(deep2_o32, bO32_DEEP, marker='3', color='r',
               label='DEEP2 Individual Spectra')
    ax.scatter(der_O32_MACT, bO32_MACT, marker='4', color='m',
               label='MACT Individual Spectra')
    for aa in range(len(valid_ID)):
        ax.annotate(valid_ID[aa], (lO32[aa], jO32_det[aa]), fontsize='6')
    ax.set_xlabel('Observed ' + r'$log(O_{32})$ ')
    ax.set_ylabel('Estimated ' + r'$log(O_{32})$')
    plt.plot(jO32_det, jO32_det, 'k', label='One to one line')
    plt.legend()

    an_txt = r'$<\Delta_{R_{23}}>$ : %0.2f' % avgO0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % medO0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sigO0
    ax.annotate(an_txt, [0.2, 0.85], xycoords='axes fraction', va='bottom',
                ha='right', fontsize=10)

    fig.savefig(pdf_pages, format='pdf')
    log.info(f"Writing: {out_pdf}")
    pdf_pages.close()

    log.info("finished.")
