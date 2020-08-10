# Calibrating my data to other work

import numpy as np
from astropy.io import ascii as asc
from os.path import exists, join

from Zcalbase_gal.Analysis import local_analog_calibration, green_peas_calibration
from Metallicity_Stack_Commons.column_names import filename_dict, npz_filename_dict

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'


def LAC_GPC_plots(fitspath, dataset, revised= False, individual = False):
    """
    Purpose
    Call function for calculating and plotting data points based with the green_pea_calibration
    and the local_analog_calibration.

    Parameters
    fitspath -> save location of the current run
    dataset  -> indicates the type of binning being used
    revised  -> indicates that revised verification table is being used
    individual -> used if individual detections from Zcalbase_gal are used

    Outputs
    pdf_files
    """
    if revised:
        temp_table = asc.read(join(fitspath,filename_dict['bin_derived_prop_rev']))
        verification = asc.read(join(fitspath, filename_dict['bin_valid_rev']))
        pea_out_pdf = join(fitspath,dataset+'_GPC.revised.pdf')
        LAC_out_pdf = join(fitspath, dataset+'_LAC.revised.pdf')
    else:
        temp_table = asc.read(join(fitspath, filename_dict['bin_derived_prop']))
        verification = asc.read(join(fitspath,filename_dict['bin_valid']))
        pea_out_pdf = join(fitspath, dataset+'_GPC.pdf')
        LAC_out_pdf = join(fitspath, dataset+'_LAC.pdf')

    detect = verification['Detection']
    det_4363 = np.where(detect == 1)[0]
    rlimit = np.where(detect == 0.5)[0]
    print('Begin Local analog Calibration')

    # Tables of individual detections from DEEP2 and MACT samples
    derived = asc.read(fitspath_ini + 'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini + 'MACT_R23_O32_derived.tbl')
    
    # DEEP2 Derived
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_OH = derived['OH'].data
    
    # MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_MACT['OH'].data

    O32_all = temp_table['logO32_Composite']
    print(O32_all)
    R23_all = temp_table['logR23_Composite']
    com_O_log = temp_table['12+log(O/H)']  # This is the 12+log(OH) value
    ID = temp_table['bin_ID']
    
    det_O32 = O32_all[det_4363]
    print('det_O32: ', det_O32)
    det_R23 = R23_all[det_4363]
    det_OH  = com_O_log[det_4363]

    rlimit_O32 = O32_all[rlimit]
    rlimit_R23 = R23_all[rlimit]
    rlimit_OH  = com_O_log[rlimit]
    rlimit_ID  = ID[rlimit]

    label = ['Detection', 'Robust Limits', 'DEEP2', 'MACT']
    marker = ['D', r'$\uparrow$', '3', '4']

    # Individual Detections from Zcalbase_gal Analysis
    if individual:
        individual_ascii = join(fitspath, filename_dict['indv_derived_prop'])
        individual = asc.read(individual_ascii)
        out_pdf = join(fitpath, 'Individual_zcalbase_gpc.pdf')
        logR23 = individual['logR23']
        logO32 = individual['logO32']
        com_log = individual['12+log(O/H)']
        bin_ID = individual['bin_ID']
        green_peas_calibration.main(logR23,logO32, com_log, out_pdf, n_bins=6, xra=[0.3,1.15],
                                    yra=[6.5,9.10], marker=['D'], edgecolors=['face', 'face', 'none'],
                                    label=['Individual Zcalbase_gal Detection'], ID = [bin_ID], fit=False,
                                    silent=False, verbose=True)
    # For LAC
    if dataset == 'R23_Grid':
        lR23 = [det_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, der_O32, der_O32_MACT]
        OH = [det_OH, der_OH, der_OH_MACT]
        local_analog_calibration.main(lR23, lO32, OH, LAC_out_pdf, ctype=['b', 'r', 'm'], label=label,
                                      silent=False, verbose=True)
        
    if dataset == 'O32_Grid' or dataset == 'Grid':    
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH   = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
        local_analog_calibration.main(lR23, lO32, OH, LAC_out_pdf, ctype=['b', 'g', 'r', 'm'],
                                      label=['Detection', 'Non-Dectection','DEEP2', 'MACT'], silent=False, verbose=True)

    if dataset == 'Voronoi10' or dataset == 'Voronoi14' or dataset == 'Voronoi20' or dataset =='Double_Bin':
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH   = [det_OH, rlimit_OH, der_OH, der_OH_MACT]

        local_analog_calibration.main(lR23, lO32, OH, LAC_out_pdf, yra=[7.0, 9.0], ctype=['b', 'g', 'r', 'm'],
                                      label=['Detection', 'Non-Dectection', 'DEEP2', 'MACT'], silent=False, verbose=True)

    if dataset == 'n_Bins':
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH   = [det_OH, rlimit_OH, der_OH, der_OH_MACT]

        local_analog_calibration.main(lR23, lO32, OH, LAC_out_pdf, yra=[7.0, 9.0], ctype=['b', 'g', 'r', 'm'],
                                      label=['Detection', 'Non-Dectection', 'DEEP2', 'MACT'], silent=False, verbose=True)
        print('finished LAC plot') 

    # For Green Pea Calibration

    if dataset == 'R23_Grid':
        lR23 = [det_R23, der_R23, der_R23_MACT]            # [det_R23,nandet_R23,der_R23,der_R23_MACT]
        lO32 = [det_O32, der_O32, der_O32_MACT]            # [det_O32,nandet_O32,der_O32,der_O32_MACT]
        OH = [det_OH, der_OH, der_OH_MACT]             # [det_OH,nandet_OH, der_OH, der_OH_MACT]

    else:
        lR23 = [det_R23, rlimit_R23, der_R23, der_R23_MACT]
        lO32 = [det_O32, rlimit_O32, der_O32, der_O32_MACT]
        OH = [det_OH, rlimit_OH, der_OH, der_OH_MACT]
        IDs = [det_ID, rlimit_ID]

        error_npz_file = join(fitspath, npz_filename_dict['der_prop_errors'])
        if exists(error_npz_file):
            print('Error npz found  ', error_npz_file, ': Adding error bars to plot' )
            error_npz = np.load(error_npz_file)
            metal_err = error_npz['12+log(O/H)_lowhigh_error'] # log values
            green_peas_calibration.main(lR23, lO32, OH, pea_out_pdf, n_bins=6, lR23_err = [], OH_err = [metal_err],
                                        xra=[0.5, 1.1], yra=[6.5, 9.10], marker=marker, label=label, IDs=IDs,
                                        include_Rlimit=True, fit=False, silent=False, verbose=True)
            
        else:
            print('No error npz found')
            green_peas_calibration.main(lR23, lO32, OH, pea_out_pdf, n_bins=6, xra=[0.5, 1.1], yra=[6.5, 9.10],
                                        marker=marker, label=label, IDs = IDs, include_Rlimit = True,
                                        fit=False, silent=False, verbose=True)


def individual_GPC(individual_ascii, validation_table):
    """
    This function is currently repetitive of the function above for the individual detection cases.
    However, I am going to keep it here in the case that we want to plot the individual detections against
    the binned detections.

    """
    pea_out_pdf_ind = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/jiang_plot_individual.pdf'
    individual = asc.read(individual_ascii)
    logR23 = individual['logR23']
    logO32 = individual['logO32']
    com_log = individual['12+log(O/H)']


    valid = asc.read(validation_table)
    Detections = valid['Detection']
    detect = np.where((Detections == 1.0))[0]
    rlimit = np.where((Detections == 0.5))[0]
    bins = valid['bin_ID']

    ID_detect = bins[detect]
    ID_rlimit = bins[rlimit]


    lR23 = [logR23]
    lO32 = [logO32]
    OH   = [com_log]
    Id = [ID_detect]

    green_peas_calibration.main(lR23, lO32, OH, pea_out_pdf_ind, n_bins=6, xra=[0.3, 1.15], yra=[6.5, 9.10],
                                marker=['3'], label=['Individual Detection'], ID=Id, fit=False,
                                silent=False, verbose=True)
