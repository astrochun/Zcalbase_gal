# This function runs the entire process start to finish
# EW values:equival width
from . import get_det3
from os.path import join
from astropy.io import ascii as asc

from . import stackboth_mastergrid, zoom_and_gauss_general, \
    calibration_plots, name_dict
from .binning import n_bins_grid_analysis, fixed_grid_analysis, \
    single_grid_o32, single_grid_r23
from .plotting import more_plots, line_ratio_plotting, te_metal_plots
from .log_commons import LogClass

from Metallicity_Stack_Commons import dir_date, valid_table, get_user
from Metallicity_Stack_Commons.column_names import filename_dict
from Metallicity_Stack_Commons.plotting import balmer
from Metallicity_Stack_Commons.analysis import error_prop


def run_grid_r23_o32_analysis(dataset, n_split=3, y_correction=False,
                              adaptive=True, apply_dust=False, mask=True,
                              thesis=False):
    """
    Function runs the entire analysis

    :param dataset: str. Define binning method options:
                    'Grid', 'O32_Grid', 'R23_Grid', 'n_Bins'
    :param n_split: int. Number of times the log(R23) bins are split.
           Default: 3.  Only for 'Grid' case
    :param y_correction: bool. Set to use smoothed ('movingaverage_box1D')
           Version for spectrum in zoom_and_gauss_general. Default: False
    :param adaptive: bool. Set to use adaptive binning, otherwise equal
           log(R23) bins. Default: True
    :param apply_dust: bool. Apply dust correction. Default: True
                       By setting apply_dust=true, all cases are run
    :param mask: bool. Mask night sky mask in stackingboth_mastergrid.
           Default: True
    :param thesis: bool. Setting that creates documents for paper writing

    No returns
    """

    fitspath_ini = get_user()
    fitspath_currentrun = join(fitspath_ini, 'Zcalbase_gal/Current_Runs/',
                               dataset)
    fitspath = dir_date(fitspath_currentrun, year=False)

    # Define logging function
    log = LogClass(fitspath, 'run_grid_r23_o32_analysis.log').get_logger()

    log.info("starting ...")

    if dataset not in ['Grid', 'O32_Grid', 'R23_Grid', 'n_Bins']:
        log.warning("Incorrect [dataset] input")
        raise ValueError("Warning!!!! Incorrect [dataset]")

    log.info(f"Run has these parameters: y_correction={y_correction}, "
             f"adaptive={adaptive}, apply_dust={apply_dust}, mask={mask}")
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

    bin_pdf_file = join(fitspath, name_dict['gridpdf_suffix'])
    bin_npz_file = join(fitspath, name_dict['gridnpz_suffix'])

    if dataset == 'O32_Grid':
        single_grid_o32.single_grid_o32(bin_pdf_file, bin_npz_file,
                                        R23, O32, galinbin, log=log)
    if dataset == 'R23_Grid':
        single_grid_r23.single_grid_r23(bin_pdf_file, bin_npz_file,
                                        R23, O32, galinbin, log=log)
    if dataset == 'Grid':
        R23_bin = 0.25
        O32_bin = 0.25
        fixed_grid_analysis.making_grid(fitspath, bin_pdf_file, bin_npz_file,
                                        R23, O32, det3, R23_bin, O32_bin,
                                        log=log)

    if dataset == 'n_Bins':
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_file,
                                            bin_npz_file, n_split,
                                            individual_ID, R23, O32,
                                            SNR3, data3, galinbin,
                                            thesis=thesis, log=log)

    log.info("made npz, pdf files, testmastergrid "
             "(need to find if this is used anywhere)")
    log.info("finished Binning_and_Graphing_MasterGrid")

    # Stackboth_MasterGrid
    # Option to Change: Bin size
    # Option to Change: Masking the night sky emission lines
    stackboth_mastergrid.master_stacking(fitspath, fitspath_ini, dataset,
                                         bin_npz_file, mask=mask, log=log)

    # Outfile and pdf both use name
    log.info(f"finished with stacking, Stacking_Masked_MasterGrid "
             f"pdf and fits files created")

    # Zoom_and_gauss_general
    # used to be 'n_Bins_2d_binning_datadet3.tbl'
    indv_bin_info = join(fitspath, filename_dict['indv_bin_info'])

    # Option to change: Constants used as initial guesses for gaussian fit

    zoom_and_gauss_general.zm_general(dataset, fitspath, y_correction, 
                                      thesis=thesis, log=log)

    log.info(f"finished gaussian fitting: {fitspath}Zoomed_Gauss_*" +
             " pdfs and fits created")
    log.info("combine_flux_table created")
    # combine_flux_ascii = join(fitspath, filename_dict['bin_fit'])

    # FIX THIS CODE: line_ratio_plotting
    # I need to go back through and figure out
    # what is the average and what is the composite
    # line_ratio_plotting.Plotting_Data1(fitspath,
    # dataset, combine_flux_ascii, binning_avg_asc)

    # Verification Table
    # If dataset is n_bins the validation table does not have to be revised
    if dataset == 'n_Bins':
        valid_table.make_validation_table(fitspath, vmin_4363SN=3.5,
                                          vmin_5007SN=95, vmax_4363sig=0.75,
                                          rlmin_4363SN=0, rlmax_4363sig=1.4,
                                          rlmin_5007SN=140)
    else:
        valid_table.make_validation_table(fitspath)

    valid_table.compare_to_by_eye(fitspath, dataset)

    # Calculating R, Temperature, Metallicity,
    # Dust Attenuation, and Errors using MSC
    # Run raw data derived properties calculations
    # (option to apply dust correction)
    error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                   apply_dust=False, revised=True, log=log)
    if apply_dust:
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=True, revised=True, log=log)

    # Run Monte Carlo randomization calculations
    # (option to apply dust correction)
    error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                   apply_dust=False, revised=True, log=log)
    if apply_dust:
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=True, revised=True, log=log)
        # use_revised here indicates weather emission.fit is .MC or not
        # emission.MC.fit is created before this so use_revised=True
        balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits',
                           use_revised=True, log=log)

    # Individual Detections
    # composite_indv_detect.main(fitspath,
    # dataset= '', revised = False, det3=True)
    # print('Individual Detections Complete')

    log.info("finished.")


def run_grid_plots(fitspath_ini, dataset, raw=False, apply_dust=False,
                   revised=False, individual=False):
    """
    This function will be run to make plots of the final results

    :param fitspath_ini: str. Path to the location of whole project
    :param dataset: str. Define binning method options:
                    'Grid', 'O32_Grid', 'R23_Grid', 'n_Bins'
    :param raw: bool. to do a simple calculation, no randomization.
                Default: False
    :param apply_dust: bool. Apply dust correction. Default: False
    :param revised: bool. indicates whether to use revised bin properties
                    (e.g., *.revised.tbl files) Default: False
    :param individual: bool. indicates whether individual galaxies are included
                       Default: False

    No returns
    """

    fitspath_currentrun = join(fitspath_ini, 'Zcalbase_gal/Current_Runs',
                               dataset)
    fitspath = dir_date(fitspath_currentrun, year=False)

    # Define logging function
    log = LogClass(fitspath, 'run_grid_plots.log').get_logger()

    log.info("starting ...")

    # Te_metal Plots
    te_metal_plots.plotting_te_metal(fitspath, fitspath_ini, raw=raw,
                                     apply_dust=apply_dust,
                                     revised=revised, individual=individual)

    # Calibration Plots
    calibration_plots.lac_gpc_plots(fitspath, fitspath_ini, dataset, raw=raw,
                                    apply_dust=apply_dust, revised=revised,
                                    individual=individual, log=log)

    '''

    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    
    temp_m_gascii = join(fitspath, 'bin_derived_properties.tbl')
    
    

    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_m_gascii, 
                           verification_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_m_gascii, 
                           verification_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_gascii, 
                                verification_table)
    more_plots.hist_for_bin(dataset, asc_table2)
    '''

    log.info("finished.")


# Below function will run the individual functions in
# the codes above that produce graphs
# Enter a keyword for want to indictate what function you want to run
# This will ease the reproduction process
# CHECK: function defaults to put new graphs in fitspath.
# Make sure you don't over write something you need
def run_individual_functions(fitspath, want, dataset='n_Bins', n_split=3,
                             adaptive=True, y_correction=False,
                             apply_dust=False, mask=True):
    """
    To run individual codes to test changes or edits plots

    :param fitspath: str. Path to the location of files saved for each run
    :param want: str. Keyword to determine what part of the process
                      needs to be run
          Options are: "binning_and_graphing", "stack_mastergrid",
                       "zoom", "R_cal_temp"
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

    No returns
    """

    log = LogClass(fitspath, 'run_individual_functions.log').get_logger()

    log.info("starting ...")

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
        bin_pdf_file = join(fitspath, name_dict['gridpdf_suffix'])
        bin_npz_file = join(fitspath, name_dict['gridnpz_suffix'])
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_file,
                                            bin_npz_file, n_split,
                                            individual_ID, R23, O32, SNR3,
                                            data3, galinbin, log=log)
        # Starting Stacking
        stackboth_mastergrid.master_stacking(fitspath, fitspath_ini, dataset,
                                             bin_npz_file, mask=mask, log=log)

    if want == 'zoom':
        zoom_and_gauss_general.zm_general(dataset, fitspath,
                                          y_correction=y_correction,
                                          log=log)

    if want == 'R_cal_temp':
        # Calculating R, Temperature, Metallicity, Dust Attenuation,
        # and Errors using MSC
        if apply_dust:
            balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits',
                               use_revised=False, log=log)

        # Run raw data derived properties calculations
        # (option to apply dust correction)
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=False, revised=False,
                                       log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                       apply_dust=False, revised=True, log=log)
        if apply_dust:
            error_prop.fluxes_derived_prop(fitspath, raw=True,
                                           binned_data=True,
                                           apply_dust=True, revised=False,
                                           log=log)
            error_prop.fluxes_derived_prop(fitspath, raw=True,
                                           binned_data=True,
                                           apply_dust=True, revised=True,
                                           log=log)

        # Run Monte Carlo randomization calculations
        # (option to apply dust correction)
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=False, revised=False,
                                       log=log)
        error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                       apply_dust=False, revised=True, log=log)
        if apply_dust:
            error_prop.fluxes_derived_prop(fitspath, raw=False,
                                           binned_data=True,
                                           apply_dust=True, revised=False,
                                           log=log)
            error_prop.fluxes_derived_prop(fitspath, raw=False,
                                           binned_data=True,
                                           apply_dust=True, revised=True,
                                           log=log)

    log.info("finished.")
