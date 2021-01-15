# This function runs the entire process start to finish
# EW values:equival width

from . import name_dict
import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack
from astropy.table import Table
from os.path import join


from . import stackboth_mastergrid, zoom_and_gauss_general, \
    calibration_plots, name_dict
from .binning import n_bins_grid_analysis, fixed_grid_analysis, \
    single_grid_o32, single_grid_r23
from .plotting import more_plots, line_ratio_plotting, te_metal_plots
from .logging import LogClass, log_stdout

from Metallicity_Stack_Commons import exclude_outliers, dir_date, lambda0, \
    valid_table, get_user
from Metallicity_Stack_Commons.column_names import filename_dict
from Metallicity_Stack_Commons.plotting import balmer
from Metallicity_Stack_Commons.analysis import error_prop


def get_det3(fitspath, fitspath_ini, log=None):
    """
    Purpose
    ----------
    Function is called to collect data for binning from DEEP2 data files

    Parameters
    ----------
    fitspath -> path to the location of files saved for each run
    fitspath_ini -> path to the location of entire project

    Returns
    -------
    Data for run
    Creates "individual_properties.tbl"
    """
    if log is None:
        log = log_stdout()

    log.info("starting ...")

    for ii in range(1, 5):
        file1 = join(fitspath_ini,
                     f"DEEP2_Commons/Catalogs/DEEP2_Field{ii}_all_line_fit.fits")
        log.info(f"Reading: {file1}")
        data = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

    objno = data0['OBJNO']

    # Excluding Outliers
    exclude_flag = exclude_outliers(objno)
    log.info(f"exclude flag: ", np.where(exclude_flag == 1))

    O2_ini = data0['OII_FLUX_MOD'].data
    O3_ini = 1.33*data0['OIIIR_FLUX_MOD'].data
    O4959_ini = data0['OIIIB_FLUX_MOD'].data
    O5007_ini = data0['OIIIR_FLUX_MOD'].data
    Hgamma_ini = data0['HG_FLUX_MOD'].data
    O4363_ini = data0['OIIIA_FLUX_MOD'].data
    Hdelta_ini = data0['HD_FLUX_MOD'].data
    
    Hb_ini = data0['HB_FLUX_MOD'].data
    R23_ini = (O2_ini+O3_ini)/Hb_ini
    O32_ini = O3_ini/O2_ini

    lR23_ini = np.log10(R23_ini)
    lO32_ini = np.log10(O32_ini)

    SNR2_ini = data0['OII_SNR'].data
    SNR3_ini = data0['OIIIR_SNR'].data
    SNRH_ini = data0['HB_SNR'].data
    SNRHg_ini = data0['HG_SNR'].data
    SNR4363_ini = data0['OIIIA_SNR'].data

    log.info(f"O2 len: {len(O2_ini)}")

    #################################################################################
    # SNR code: This rules out major outliers by only using specified data
    # May limit the logR23 value further to 1.2, check the velocity dispersions of the high R23 spectra
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini > 0) &
                    (exclude_flag == 0) & (lR23_ini < 1.4))[0]

    # Organize the R23_032 data
    data3 = data0[det3]

    R23 = R23_ini[det3]
    O32 = O32_ini[det3]
    lR23 = lR23_ini[det3]
    lO32 = lO32_ini[det3]
    
    Hb = Hb_ini[det3]
    O2 = O2_ini[det3]
    O3 = O3_ini[det3]
    Hgamma = Hgamma_ini[det3]
    Hdelta = Hdelta_ini[det3]
    O4363 = O4363_ini[det3]
    O4959 = O4959_ini[det3]
    O5007 = O5007_ini[det3]
    SNR2 = SNR2_ini[det3]
    SNR3 = SNR3_ini[det3]
    SNRH = SNRH_ini[det3]
    SNRHG = SNRHg_ini[det3]
    SNR4363 = SNR4363_ini[det3]
    individual_names = objno[det3]

    table_dict = {
        'ID': individual_names,
        'logR23': lR23,
        'logO32': lO32,
        'OII_3727_Flux_Gaussian': O2,
        'O3_Flux_Gaussian': O3,
        'HGAMMA_Flux_Gaussian': Hgamma,
        'HDELTA_Flux_Gaussian': Hdelta,
        'OIII_4363_Flux_Gaussian': O4363,
        'OIII_4958_Flux_Gaussian': O4959,
        'OIII_5007_Flux_Gaussian': O5007,
        'HBETA_Flux_Gaussian': Hb,
        'O2_S/N': SNR2,
        'O3_S/N': SNR3,
        'RH_S/N': SNRH,
        'HGAMMA_S/N': SNRHG,
        'O4363_S/N': SNR4363
    }

    tab1 = Table(table_dict)

    # We can create two different kinds of tables here of the R23_032 data (det3)
    # used to be get_det3_table.tbl
    outfile = join(fitspath, filename_dict['indv_prop'])
    log.info(f"Writing: {outfile}")
    asc.write(tab1, outfile, format='fixed_width_two_line')
    # tab1.write(fitspath_ini+'get_det3_table.fit', format = 'fits', overwrite = True)

    log.info("finished.")

    return individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, data3


def run_grid_r23_o32_analysis(dataset, n_split=3, y_correction=False,
                              adaptive=True, apply_dust=False, mask=True):
    """
    Purpose
      Function runs the entire analysis

    :param dataset: str. Define binning method options:
                    'Grid', 'O32_Grid', 'R23_Grid', 'n_Bins'
    :param n_split: int. Number of times the log(R23) bins are split.
           Default: 3.  Only for 'Grid' case
    :param y_correction: bool. Set to use smoothed ('movingaverage_box1D')
           Version for spectrum in zoom_and_gauss_general. Default: False
    :param adaptive: bool. Set to use adaptive binning, otherwise equal
           log(R23) bins. Default: True
    :param apply_dust: bool. Apply dust correction. Default: False
    :param mask: bool. Mask night sky mask in stackingboth_mastergrid.
           Default: True
    """

    fitspath_ini = get_user()
    fitspath_currentrun = join(fitspath_ini, 'Zcalbase_gal/Current_Runs/')
    fitspath = dir_date(fitspath_currentrun, year=False)

    # Define logging function
    log = LogClass(fitspath, 'run_grid_r23_o32_analysis.log').get_logger()

    if dataset not in ['Grid', 'O32_Grid', 'R23_Grid', 'n_Bins']:
        log.warning("Incorrect [dataset] input")
        raise ValueError("Warning!!!! Incorrect [dataset]")

    log.info(f"fitspath_ini = {fitspath_ini}")
    log.info(f"fitspath = {fitspath}")

    individual_ID, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, \
        data3 = get_det3(fitspath, fitspath_ini, log=log)

    log.info(f"length R23: {len(R23)}")

    # Each bin will be split in half
    # Must sum to 2799 for Zcalbase_gal Analysis
    if adaptive:
        galinbin = [458, 450, 400, 300, 300, 275, 250, 200, 176]
    else:
        galinbin = [400, 400, 400, 400, 400, 400, 409]
    log.info(f"# of Gal in Bin: {galinbin}")

    bin_pdf_pages = join(fitspath, dataset + name_dict['gridpdf_suffix'])
    bin_outfile = join(fitspath, dataset + name_dict['gridnpz_suffix'])

    if dataset == 'O32_Grid':
        single_grid_o32.single_grid_o32(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, galinbin, log=log)
    if dataset == 'R23_Grid':
        single_grid_r23.single_grid_r23(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, galinbin, log=log)
    if dataset == 'Grid':
        R23_bin = 0.25
        O32_bin = 0.25
        fixed_grid_analysis.making_grid(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, det3, R23_bin, O32_bin, log=log)

    if dataset == 'n_Bins':
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_pages, bin_outfile,
                                            n_split, individual_ID, R23, O32,
                                            SNR3, data3, galinbin, log=log)

    log.info("made npz, pdf files, testmastergrid (need to find if this is used anywhere)")
    log.info("finished Binning_and_Graphing_MasterGrid")

    # Stackboth_MasterGrid
    # Option to Change: Bin size
    # Option to Change: Masking the night sky emission lines

    if mask:
        stack_name = join(dataset, name_dict['Stackname'])
    else:
        stack_name = join(dataset, name_dict['Stackname_nomask'])

    stackboth_mastergrid.master_stacking(fitspath, fitspath_ini, dataset,
                                         bin_outfile, stack_name, mask=mask,
                                         log=log)

    # Outfile and pdf both use name
    log.info(f"finished with stacking, {stack_name} pdf and fits files created")

    # Zoom_and_gauss_general
    outfile_grid = join(fitspath, filename_dict['comp_spec'])
    log.info(f"outfile_grid : {outfile_grid}")
    log.info(f"Reading: {outfile_grid}")
    stack2D, header = fits.getdata(outfile_grid, header=True)
    wave = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
    dispersion = header['CDELT1']
    binning_avg_asc = join(fitspath, filename_dict['bin_info'])

    # used to be 'n_Bins_2d_binning_datadet3.tbl'
    indv_bin_info = join(fitspath, filename_dict['indv_bin_info'])

    lineflag = np.zeros(len(wave))
    for ii in lambda0:   
        idx = np.where(np.absolute(wave - ii) <= 5)[0]
        if len(idx) > 0:
            lineflag[idx] = 1
    # Option to change: Constants used as initial guesses for gaussian fit

    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave, lineflag,
                                      dispersion, y_correction=y_correction,
                                      tab=binning_avg_asc, log=log)

    log.info(f"finished gaussian fitting: {fitspath}_{dataset}_Zoomed_Gauss_*" +
             " pdfs and fits created")
    log.info("combine_flux_table created")
    # combine_flux_ascii = join(fitspath, filename_dict['bin_fit'])

    # FIX THIS CODE: line_ratio_plotting
    # I need to go back through and figure out what is the average and what is the composite
    # line_ratio_plotting.Plotting_Data1(fitspath, dataset, combine_flux_ascii, binning_avg_asc)

    # Verification Table
    valid_table.make_validation_table(fitspath)
    verification_table = join(fitspath, filename_dict['bin_valid'])

    valid_table.compare_to_by_eye(fitspath, dataset)
    verification_table_revised = join(fitspath, filename_dict['bin_valid_rev'])

    # Calculating R, Temperature, Metallicity, Dust Attenuation, and Errors using MSC
    if apply_dust:
        balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits',
                           use_revised=False, log=log)

    # Run raw data derived properties calculations (option to apply dust correction)
    error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                   apply_dust=False, revised=False, log=log)
    error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                   apply_dust=False, revised=True, log=log)
    if apply_dust:
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=True, revised=False, log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=True, revised=True, log=log)

    # Run Monte Carlo randomization calculations (option to apply dust correction)
    error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                   apply_dust=False, revised=False, log=log)
    error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                   apply_dust=False, revised=True, log=log)
    if apply_dust:
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=True, revised=False, log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=True, revised=True, log=log)

    # Individual Detections
    # composite_indv_detect.main(fitspath, dataset= '', revised = False, det3=True)
    # print('Individual Detections Complete')


def run_grid_plots(fitspath, fitspath_ini, dataset, raw=False, apply_dust=False,
                   revised=False, individual=False):
    """
    Purpose
    --------
    This function will be run to make plots of the final results

    Parameters
    ----------

    """
    # Te_metal Plots
    te_metal_plots.plotting_te_metal(fitspath, fitspath_ini, raw=raw, apply_dust=apply_dust,
                                     revised=revised, individual=individual)

    # Calibration Plots
    calibration_plots.lac_gpc_plots(fitspath, fitspath_ini, dataset, raw=raw,
                                    apply_dust=apply_dust, revised=revised,
                                    individual=individual, log=log)

    '''

    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    #asc_table_det3 = asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    
    temp_m_gascii = join(fitspath, 'bin_derived_properties.tbl')
    
    

    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_m_gascii, verification_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_m_gascii, verification_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_gascii, verification_table)
    more_plots.hist_for_bin(dataset, asc_table2)
    '''

    log.debug("finished.")


# Below function will run the individual functions in the codes above that produce graphs
# Enter a keyword for want to indictate what function you want to run
# This will ease the reproduction process
# CHECK: function defaults to put new graphs in fitspath. Make sure you don't over write something you need
def run_individual_functions(fitspath, want, dataset='n_Bins', n_split=3,
                             adaptive=True, y_correction=False,
                             apply_dust=False, mask=True):
    """
    Purpose
      To run individual codes to test changes or edits plots

    Parameters
    ----------
    :param fitspath: str. Path to the location of files saved for each run
    :param want: str. Keyword to determine what part of the process needs to be run
          Options are: "binning_and_graphing", "stack_mastergrid", "zoom", "R_cal_temp"
    :param dataset: str. Define binning method options:
                    'Grid', 'O32_Grid', 'R23_Grid', 'n_Bins'
    :param n_split: int. Number of times the log(R23) bins are split.
           Default: 3.  Only for 'Grid' case
    :param y_correction: bool. Set to use smoothed ('movingaverage_box1D')
           Version for spectrum in zoom_and_gauss_general. Default: False
    :param adaptive: bool. Set to use adaptive binning, otherwise equal
           log(R23) bins. Default: True
    :param apply_dust: bool. Apply dust correction. Default: False
    :param mask: bool. Mask night sky mask in stackingboth_mastergrid.
           Default: True
    """

    log = LogClass(fitspath, 'run_individual_functions.log').get_logger()

    fitspath_ini = get_user()
    if want == 'binning_and_stacking':
        individual_ID, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, \
            data3 = get_det3(fitspath, fitspath_ini, log=log)
        # Each bin will be split in half
        # Must sum to 2799
        if adaptive:
            galinbin = [458, 450, 400, 300, 300, 275, 250, 200, 176]
        else:
            galinbin = [400, 400, 400, 400, 400, 400, 409]
        bin_pdf_pages = join(fitspath, name_dict['gridpdf_suffix'])
        bin_outfile = join(fitspath, name_dict['gridnpz_suffix'])
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_pages,
                                            bin_outfile, n_split,
                                            individual_ID, R23, O32, SNR3,
                                            data3, galinbin, log=log)
        # Starting Stacking
        if mask:
            stack_name = dataset + name_dict['Stackname']
        else:
            stack_name = dataset + name_dict['Stackname_nomask']
        stackboth_mastergrid.master_stacking(fitspath, fitspath_ini, dataset,
                                             bin_outfile, stack_name,
                                             mask=mask, log=log)

    if want == 'zoom':
        outfile_grid = join(fitspath, filename_dict['comp_spec'])
        log.info(f"Reading: {outfile_grid}")
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        dispersion = header['CDELT1']
        binning_avg_asc = join(fitspath, filename_dict['bin_info'])

        lineflag = np.zeros(len(wave))
        for ii in lambda0:   
            idx = np.where(np.absolute(wave - ii) <= 5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave,
                                          lineflag, dispersion,
                                          y_correction=y_correction,
                                          tab=binning_avg_asc, log=log)

    if want == 'R_cal_temp':
        # Calculating R, Temperature, Metallicity, Dust Attenuation, and Errors using MSC
        if apply_dust:
            balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits',
                               use_revised=False, log=log)

        # Run raw data derived properties calculations (option to apply dust correction)
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=False, revised=False, log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=False, revised=True, log=log)
        if apply_dust:
            error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                           apply_dust=True, revised=False, log=log)
            error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                           apply_dust=True, revised=True, log=log)

        # Run Monte Carlo randomization calculations (option to apply dust correction)
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=False, revised=False, log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=False, revised=True, log=log)
        if apply_dust:
            error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                           apply_dust=True, revised=False, log=log)
            error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                           apply_dust=True, revised=True, log=log)

    log.debug("finished.")
