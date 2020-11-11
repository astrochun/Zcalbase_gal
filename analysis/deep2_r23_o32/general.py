# This function runs the entire process start to finish
# This weekend combine the grid and voronoi if statements, voronoi20, log plots for met
# EW values:equival width

from . import name_dict
import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack
from astropy.table import Table
from os.path import join


from Zcalbase_gal.analysis.deep2_r23_o32 import stackboth_mastergrid, \
    zoom_and_gauss_general, hstack_tables, calibration_plots, name_dict
from Zcalbase_gal.analysis.deep2_r23_o32.binning import n_bins_grid_analysis, fixed_grid_analysis, \
    single_grid_o32, single_grid_r23
from Zcalbase_gal.analysis.deep2_r23_o32.plotting import more_plots, line_ratio_plotting, te_metal_plots
from Metallicity_Stack_Commons import exclude_outliers, dir_date, lambda0, \
    valid_table, get_user
from Metallicity_Stack_Commons.column_names import filename_dict
from Metallicity_Stack_Commons.plotting import balmer
from Metallicity_Stack_Commons.analysis import error_prop


def get_det3(fitspath, fitspath_ini):
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
    for ii in range(1, 5):
        file1 = fitspath_ini+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

    objno = data0['OBJNO']

    # Excluding Outliers
    exclude_flag = exclude_outliers(objno)
    print("exclude flag: ", np.where(exclude_flag == 1))

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

    print('O2 len:', len(O2_ini))

    #################################################################################
    # SNR code: This rules out major outliers by only using specified data
    # May limit the logR23 value further to 1.2, check the velocity dispersions of the high R23 spectra
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini > 0) & (exclude_flag == 0) & (lR23_ini < 1.4))[0]

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
    asc.write(tab1, join(fitspath, filename_dict['indv_prop']),
              format='fixed_width_two_line')
    # tab1.write(fitspath_ini+'get_det3_table.fit', format = 'fits', overwrite = True)

    return individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3


def run_grid_r23_o32_analysis(dataset, n_split=3, y_correction=False, adaptive=True,
                              apply_dust=False, mask=True):
    """
    Purpose
    ----------
    Function runs the entire analysis

    Parameters
    ----------
    :param dataset: str. Define binning method options:
                    'Grid', 'O32_Grid', 'R23_Grid', 'n_Bins'
    :param n_split: int. Number of times the log(R23) bins are split.
           Default: 3.  Only for 'Grid' case
    :param y_correction: bool. Set to use smoothed (movingaverage_box1D)
           version for spectrum in zoom_and_gauss_general. Default: False
    :param adaptive: bool. Set to use adaptive binning, otherwise equal
           log(R23) bins. Default: True
    :param apply_dust: bool. Apply dust correction. Default: False
    :param mask: bool. Mask night sky mask in stackingboth_mastergrid.
           Default: True
    """

    if dataset not in ['Grid', 'O32_Grid', 'R23_Grid', 'n_Bins']:
        print("Warning!!!! Incorrect [dataset]")
        raise ValueError("Warning!!!! Incorrect [dataset]")

    fitspath_ini = get_user()
    fitspath = dir_date(fitspath_ini, year=False)
    print('fitspath_ini =  ', fitspath_ini)
    print('fitspath =  ', fitspath)

    individual_ID, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, \
        data3 = get_det3(fitspath, fitspath_ini)

    print("length R23: ", len(R23))

    # Each bin will be split in half
    # Must sum to 2799 for Zcalbase_gal Analysis
    if adaptive:
        galinbin = [458, 450, 400, 300, 300, 275, 250, 200, 176]
    else:
        galinbin = [400, 400, 400, 400, 400, 400, 409]
    print('# of Gal in Bin:', galinbin)

    bin_pdf_pages = join(fitspath, dataset + name_dict['gridpdf_suffix'])
    bin_outfile = join(fitspath, dataset + name_dict['gridnpz_suffix'])

    if dataset == 'O32_Grid':
        single_grid_o32.single_grid_o32(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH,
                                        det3, data3, galinbin, adaptive)
    if dataset == 'R23_Grid':
        single_grid_r23.single_grid_r23(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH,
                                        det3, data3, galinbin)
    if dataset == 'Grid':
        R23_bin = 0.25
        O32_bin = 0.25
        fixed_grid_analysis.making_Grid(fitspath, bin_pdf_pages, bin_outfile,
                                        R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH,
                                        det3, data3, R23_bin, O32_bin)

    if dataset == 'n_Bins':
        n_bins_grid_analysis.n_times_binned(fitspath, bin_pdf_pages, bin_outfile,
                                            n_split, individual_ID, R23, O32,
                                            SNR3, data3, galinbin)

    print('made npz, pdf files , testmastergrid (need to find if this is used anywhere)')
    print('finished Binning_and_Graphing_MasterGrid')

    # Stackboth_MasterGrid
    # Option to Change: Bin size
    # Option to Change: Masking the night sky emission lines

    if mask:
        stack_name = dataset + name_dict['Stackname']
    else:
        stack_name = dataset + name_dict['Stackname_nomask']

    stackboth_mastergrid.master_stacking(fitspath, fitspath_ini, dataset,
                                         bin_outfile, stack_name, mask=mask)

    # Outfile and pdf both use name
    print('finished with stacking,' + stack_name + 'pdf and fits files created')

    # Zoom_and_gauss_general
    outfile_grid = join(fitspath, filename_dict['comp_spec'])
    print(outfile_grid)
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
                                      tab=binning_avg_asc)

    print('finished gaussian fitting:,' + fitspath + '_'+dataset+'_Zoomed_Gauss_* pdfs and fits created')
    print('combine_flux_table created')
    # combine_flux_ascii = join(fitspath, filename_dict['bin_fit'])

    # FIX THIS CODE: line_ratio_plotting
    # I need to go back through and figure out what is the average and what is the composite
    # line_ratio_plotting.Plotting_Data1(fitspath, dataset, combine_flux_ascii, binning_avg_asc)

    # Calculating R, Temperature, Metallicity, Dust Attenuation, and Errors using MSC
    if apply_dust:
        balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits', use_revised=False)

    # Generate baseline table (raw=True)
    error_prop.fluxes_derived_prop(fitspath, raw=True, binned_data=True,
                                   apply_dust=apply_dust, revised=True)

    # Generate MC randomized table (raw=False)
    error_prop.fluxes_derived_prop(fitspath, raw=False, binned_data=True,
                                   apply_dust=apply_dust, revised=True)

    '''
    ttab = asc.read(temp)
    trtab = asc.read(temp_revised)
    
    metallicity = ttab['12+log(O/H)'].data
    metallicity_r = trtab['12+log(O/H)'].data

    print('Metallicity' , metallicity)
    print('##################################################################')
    print('Dust Attenuated Metallicity', metallicity_r)'''

    # Individual Detections
    # composite_indv_detect.main(fitspath, dataset= '', revised = False, det3=True)
    # print('Individual Detections Complete')

    # Te_metal Plots
    # te_metal_plots.plotting_te_metal(fitspath, revised=False)
    # te_metal_plots.plotting_te_metal(fitspath,   revised=True)

    # Calibration Plots
    # calibration_plots.LAC_GPC_plots(fitspath, dataset, revised= False)
    calibration_plots.lac_gpc_plots(fitspath, fitspath_ini, dataset,
                                    revised=True, individual=False)


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


# Below function will run the individual functions in the codes above that produce graphs
# Enter a keyword for want to indictate what function you want to run
# This will ease the reproduction process
# CHECK: function defaults to put new graphs in fitspath. Make sure you don't over write something you need
def run_individual_functions(fitspath, want, adaptive, y_correction=False,
                             apply_dust=False):
    """
    Purpose
    ----------
    To run individual codes to test changes or edits plots

    Parameters
    ----------
    fitspath -> path to the location of files saved for each run
    want -> keyword to determine what part of the process needs to be run
         -> Keywords: binning_and_graphing, stack_mastergrid, zoom, R_cal_temp, line_ratio_plotting
    adaptive  -> determines if the R23 bins have equal or different number of spectra in them in binning method
    y_correction -> determines if the smoothed (movingaverage_box1D) version of y is used in zoom_and_gauss_general.py
    apply_dust -> determines if dust attenuation corrections are applied

    """

    dataset = 'n_Bins'

    if want == 'binning_and_graphing':
        R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3(fitspath)
        # Each bin will be split in half
        # Must sum to 2799
        if adaptive:
            galinbin = [458, 450, 400, 300, 300, 275, 250, 200, 176]
        else:
            galinbin = [400, 400, 400, 400, 400, 400, 409]
        pdf_pages = join(fitspath, 'n_Bins_grid.pdf')
        grid_data_file = join(fitspath, 'n_Bins_grid.npz')
        n_bins_grid_analysis.n_times_binned(fitspath, pdf_pages, grid_data_file, n_split,
                                            R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH,
                                            det3, data3, galinbin, adaptive)

    if want == 'stack_mastergrid':
        grid_data_file = join(fitspath, 'n_Bins_grid.npz')
        Stack_name = 'Stacking_Masked_MasterGrid_' + dataset + '.pdf'
        stackboth_masterGrid.run_Stacking_Master_mask(det3, data3, fitspath, fitspath_ini,
                                                      dataset, Stack_name, grid_data_file)

    if want == 'zoom':
        Stack_name = 'Stacking_Masked_MasterGrid_' + dataset + '.fits'
        outfile_grid = join(fitspath, Stack_name)
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        dispersion = header['CDELT1']
        binning_avg_asc = join(fitspath, 'bin_info.tbl')

        lineflag = np.zeros(len(wave))
        for ii in lambda0:   
            idx = np.where(np.absolute(wave - ii) <= 5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave,
                                          lineflag, dispersion,
                                          y_correction=y_correction,
                                          tab=binning_avg_asc)

    if want == 'R_cal_temp':
        combine_flux_ascii = join(fitspath, 'bin_emission_line_fit.tbl')
        temp_m_gascii = join(fitspath, 'nsplit_temperatures_metalicity.tbl')
        temp_m_gfits = join(fitspath, 'nsplit_temperatures_metalicity.fits')
        temp_m_gpdf_name = 'nsplit_Temp_Composite_Metallicity.pdf'

        if apply_dust:
            r_temp_calcul.run_function(fitspath, dataset, temp_m_gascii, temp_m_gfits,
                                       temp_m_gpdf_name, combine_flux_ascii, apply_dust=True)
        else:
            r_temp_calcul.run_function(fitspath, dataset, temp_m_gascii, temp_m_gfits,
                                       temp_m_gpdf_name, combine_flux_ascii, apply_dust=False)

    if want == 'line_ratio_plotting':
        combine_flux_ascii = join(fitspath, 'bin_emission_line_fit.tbl')
        binning_avg_asc = join(fitspath, 'bin_info.tbl')
        line_ratio_plotting.Plotting_Data1(fitspath, dataset, combine_flux_ascii, binning_avg_asc)

    if want == 'calibration_plots':
        temp_m_gascii = join(fitspath, 'nsplit_temperatures_metalicity.tbl')
        calibration_plots.LAC_GPC_plots(fitspath, dataset, temp_m_gascii)

    print(want, 'is done')
