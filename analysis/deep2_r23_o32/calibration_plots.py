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

    log.info("starting ...")

    suffix = ''
    if not revised:
        suffix += '.valid1'

    if not raw:
        suffix += '.MC'

    if apply_dust:
        suffix += '.dustcorr'

    gpc_pdf_file = join(fitspath, f"GPC{suffix}.pdf")
    lac_pdf_file = join(fitspath, f"LAC{suffix}.pdf")
    lac_diff_R23_pdf_file = join(fitspath, f"LAC_R23_diff{suffix}.pdf")
    lac_diff_O32_pdf_file = join(fitspath, f"LAC_O32_diff{suffix}.pdf")

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
    log.debug(lO32_all)
    com_O_log = bin_derived_prop_tab['12+log(O/H)']  # This is the 12+log(OH) value

    det_ID = ID[det_4363]
    det_lR23 = lR23_all[det_4363]
    det_lO32 = lO32_all[det_4363]
    log.debug(f"det_lO32: {det_lO32}")
    det_OH = com_O_log[det_4363]

    rlimit_ID = ID[rlimit]
    rlimit_lR23 = lR23_all[rlimit]
    rlimit_lO32 = lO32_all[rlimit]
    rlimit_OH = com_O_log[rlimit]

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
        lR23 = [det_lR23, DEEP2_lR23, MACT_lR23]
        lO32 = [det_lO32, DEEP2_lO32, MACT_lO32]
        OH = [det_OH, DEEP2_OH, MACT_OH]
        c_var = ['b', 'r', 'm']
        label = ['Detection', 'DEEP2', 'MACT']
        marker = ['D', '3', '4']
        
    if dataset in ['O32_Grid', 'Grid']:
        lR23 = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32 = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']

    if dataset in ['Voronoi10', 'Voronoi14', 'Voronoi20', 'Double_Bin']:
        lR23 = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32 = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']

    if dataset == 'n_Bins':
        lR23 = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
        lO32 = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
        OH = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
        c_var = ['b', 'g', 'r', 'm']
        label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']
        IDs = [det_ID, rlimit_ID, DEEP2_id, MACT_ID]

    '''local_analog_calibration.main(lR23, lO32, OH, lac_pdf_file,
                                  lac_diff_R23_pdf_file,
                                  ID=IDs, yra=[7.0, 9.0],
                                  ctype=c_var, label=label, marker=marker,
                                  log=log)'''

    log.info('finished LAC plot')

    # For Green Pea Calibration
    lR23 = [det_lR23, rlimit_lR23, DEEP2_lR23, MACT_lR23]
    lO32 = [det_lO32, rlimit_lO32, DEEP2_lO32, MACT_lO32]
    OH = [det_OH, rlimit_OH, DEEP2_OH, MACT_OH]
    IDs = [det_ID, rlimit_ID]

    if dataset == 'R23_Grid':
        lR23 = [det_lR23, DEEP2_lR23, MACT_lR23]
        lO32 = [det_lO32, DEEP2_lO32, MACT_lO32]
        OH = [det_OH, DEEP2_OH, MACT_OH]
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
                                        xra=[0.5, 1.1], yra=[7.0, 8.8],
                                        marker=marker, edgecolors=edgecolor,
                                        alpha=alpha, label=label, IDs=IDs,
                                        include_Rlimit=True, fit=False,
                                        log=log)
        else:
            log.info('No error npz found')
            green_peas_calibration.main(lR23, lO32, OH, gpc_pdf_file, n_bins=6,
                                        xra=[0.5, 1.1], yra=[7.0, 8.8],
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
