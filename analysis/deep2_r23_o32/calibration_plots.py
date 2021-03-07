import numpy as np
from astropy.io import ascii as asc
from os.path import exists, join

from .. import local_analog_calibration, green_peas_calibration
from Metallicity_Stack_Commons.column_names import filename_dict, \
    npz_filename_dict

from .log_commons import log_stdout


def lac_gpc_plots(fitspath, fitspath_ini, dataset, raw=False,
                  apply_dust=False, revised=False, individual=False,
                  log=None):
    """
    Call function for calculating and plotting data points based with
    the green_pea_calibration and the local_analog_calibration.

    :param fitspath: str. save location of the current run
    :param fitspath_ini: str.
    :param dataset: str. indicates the type of binning being used
    :param apply_dust: bool. Apply dust correction. Default: False
    :param revised: bool. indicates that revised verification table
                    is being used
    :param individual: bool.used if individual detections
                    from Zcalbase_gal are used
    :param log: LogClass or logging object

    PDF Files:
        fitspath + GPC{suffix}.pdf
        fitspath + LAC{suffix}.pdf
        fitspath + GPC{suffix}.diff.pdf
        fitspath + LAC{suffix}.diff.pdf
    No returns
    """

    if log is None:
        log = log_stdout()

    # log.info("starting ...")

    suffix = ''
    if not revised:
        suffix += '.valid1'

    if not raw:
        suffix += '.MC'

    if apply_dust:
        suffix += '.dustcorr'

    gpc_pdf_file = join(fitspath, f"GPC{suffix}.pdf")
    lac_pdf_file = join(fitspath, f"LAC{suffix}.pdf")
    R23_pdf_diff = join(fitspath, f"LAC{suffix}_diff.pdf")
    O32_pdf_diff = join(fitspath, f"LAC{suffix}_diff.pdf")

    # Validation Table Call
    if revised:
        bin_valid_file = join(fitspath, filename_dict['bin_valid_rev'])
    else:
        bin_valid_file = join(fitspath, filename_dict['bin_valid'])
    log.info(f"Reading: {bin_valid_file}")
    bin_valid_tab = asc.read(bin_valid_file)

    bin_derived_prop_file = join(fitspath, f"bin_derived_properties{suffix}.tbl")
    log.info(f"Reading: {bin_derived_prop_file}")
    bin_derived_prop_tab = asc.read(bin_derived_prop_file)

    detect = bin_valid_tab['Detection']
    det_4363 = np.where(detect == 1)[0]
    rlimit = np.where(detect == 0.5)[0]
    # print('det_4363: ', det_4363)
    # print('Begin Local analog Calibration')

    # Tables of individual detections from DEEP2 and MACT samples
    derived_deep_file = join(fitspath_ini, 'DEEP2_Commons/Catalogs/',
                             'DEEP2_R23_O32_derived.tbl')
    log.info(f"Reading: {derived_deep_file}")
    derived_deep_tab = asc.read(derived_deep_file)

    derived_mact_file = join(fitspath_ini, 'MACT_Commons/Catalogs/',
                             'MACT_R23_O32_derived.tbl')
    log.info(f"Reading: {derived_mact_file}")
    derived_mact_tab = asc.read(derived_mact_file)

    # DEEP2 Derived
    er_R23 = derived_deep_tab['R23'].data
    er_O32 = derived_deep_tab['O32'].data
    der_ID = derived_deep_tab['ID'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_OH = derived_deep_tab['OH'].data

    # MACT Derived
    er_R23_MACT = derived_mact_tab['R23'].data
    er_O32_MACT = derived_mact_tab['O32'].data
    der_ID_MACT = derived_mact_tab['ID'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_mact_tab['OH'].data

    O32_all = bin_derived_prop_tab['logO32_composite']
    log.debug(O32_all)
    R23_all = bin_derived_prop_tab['logR23_composite']
    com_O_log = bin_derived_prop_tab['12+log(O/H)']  # This is the 12+log(OH) value
    ID = bin_derived_prop_tab['bin_ID']
    
    det_O32 = O32_all[det_4363]
    log.debug(f"det_O32: {det_O32}")
    det_R23 = R23_all[det_4363]
    det_OH = com_O_log[det_4363]
    det_ID = ID[det_4363]

    rlimit_O32 = O32_all[rlimit]
    rlimit_R23 = R23_all[rlimit]
    rlimit_OH = com_O_log[rlimit]
    rlimit_ID = ID[rlimit]

    label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']
    marker = ['D', r'$\uparrow$', '3', '4']

    # Individual Detections from Zcalbase_gal Analysis
    if individual:
        indv_derived_prop_file = join(fitspath,
                                      filename_dict['indv_derived_prop'])
        log.info(f"Reading: {indv_derived_prop_file}")
        indv_derived_prop_tab = asc.read(indv_derived_prop_file)

        pdf_file = join(fitspath, 'Individual_zcalbase_gpc.pdf')
        logR23 = indv_derived_prop_tab['logR23']
        logO32 = indv_derived_prop_tab['logO32']
        com_log = indv_derived_prop_tab['12+log(O/H)']
        bin_ID = indv_derived_prop_tab['bin_ID']
        alpha = [1]
        green_peas_calibration.main(logR23, logO32, com_log, pdf_file, n_bins=6,
                                    xra=[0.3, 1.15], yra=[6.5, 9.10],
                                    marker=['D'],
                                    edgecolors=['face', 'face', 'none'],
                                    alpha=alpha, IDs=[bin_ID],
                                    label=['Individual Zcalbase_gal Detection'],
                                    fit=False)

    # For LAC
    if dataset == 'R23_Grid':
        lR23 = [det_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, der_O32, der_O32_MACT]
        OH = [det_OH, der_OH, der_OH_MACT]
        c_var = ['b', 'r', 'm']
        label = ['Detection', 'DEEP2', 'MACT']
        marker = ['D', '3', '4']
        
    if dataset in ['O32_Grid', 'Grid']:
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']

    if dataset in ['Voronoi10', 'Voronoi14', 'Voronoi20', 'Double_Bin']:
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']

    if dataset == 'n_Bins':
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']
        IDs = [det_ID, rlimit_ID, der_ID, der_ID_MACT]

    local_analog_calibration.main(lR23, lO32, OH, lac_pdf_file,
                                  R23_pdf_diff, O32_pdf_diff, ID=IDs,
                                  yra=[7.0, 9.0],
                                  ctype=c_var, label=label, marker=marker,
                                  log=None)

    log.info('finished LAC plot')

    # For Green Pea Calibration
    lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
    lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
    OH = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
    IDs = [det_ID, rlimit_ID]

    if dataset == 'R23_Grid':
        lR23 = [det_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, der_O32, der_O32_MACT]
        OH = [det_OH, der_OH, der_OH_MACT]
        IDs = [det_ID]

    if revised:
        alpha = np.repeat(1, len(lR23))
        edgecolor = np.repeat('face', len(lR23))
        error_npz_file = join(fitspath, npz_filename_dict['der_prop_errors'])
        if exists(error_npz_file):
            log.info(f"Error npz found {error_npz_file}: "
                     f"Adding error bars to plot")
            error_npz = np.load(error_npz_file)
            metal_err = error_npz['12+log(O/H)_error']  # log values
            green_peas_calibration.main(lR23, lO32, OH, gpc_pdf_file, n_bins=6,
                                        lR23_err=[], OH_err=[metal_err],
                                        xra=[0.5, 1.1], yra=[6.5, 9.10],
                                        marker=marker, edgecolors=edgecolor,
                                        alpha=alpha, label=label, IDs=IDs,
                                        include_Rlimit=True, fit=False,
                                        log=log)
        else:
            log.info('No error npz found')
            green_peas_calibration.main(lR23, lO32, OH, gpc_pdf_file, n_bins=6,
                                        xra=[0.5, 1.1], yra=[6.5, 9.10],
                                        marker=marker, edgecolors=edgecolor,
                                        alpha=alpha, label=label, IDs=IDs,
                                        include_Rlimit=True, fit=False,
                                        log=log)

    log.info("finished.")


def individual_gpc(indv_file, valid_file, name, log=None):
    """
    This function is currently repetitive of the function above
    for the individual detection cases.
    However, I am going to keep it here in the case that we want
    to plot the individual detections against
    the binned detections.

    :param indv_file: str. location of table with data
    :param valid_file: str. location of validaition table
    :param name: str. name of the outputted pdf file
    :param log: LogClass or logging object

    PDF Files:
        fitspath + GPC{suffix}.pdf
        fitspath + GPC{suffix}.diff.pdf

    No returns
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    log.info(f"Reading: {indv_file}")
    indv_tab = asc.read(indv_file)
    logR23 = indv_tab['logR23']
    logO32 = indv_tab['logO32']
    com_log = indv_tab['12+log(O/H)']

    log.info(f"Reading: {valid_file}")
    valid_tab = asc.read(valid_file)
    Detections = valid_tab['Detection']
    detect = np.where((Detections == 1.0))[0]
    rlimit = np.where((Detections == 0.5))[0]
    bins = valid_tab['bin_ID']

    ID_detect = bins[detect]
    ID_rlimit = bins[rlimit]

    lR23 = [logR23]
    lO32 = [logO32]
    OH = [com_log]
    Id = [ID_detect]

    green_peas_calibration.main(lR23, lO32, OH, name, n_bins=6,
                                xra=[0.3, 1.15], yra=[6.5, 9.10],
                                marker=['3'], label=['Individual Detection'],
                                IDs=Id, fit=False, log=log)
