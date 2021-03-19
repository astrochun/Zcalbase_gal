from astropy.io import fits
import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table, vstack
from os.path import join, exists

from .log_commons import log_stdout

from Metallicity_Stack_Commons import exclude_outliers
from Metallicity_Stack_Commons.column_names import filename_dict

# Dictionary of pdf and other files specific to the Zcalbase_gal project
name_dict = dict()

name_dict['gridnpz_suffix'] = 'grid.npz'
name_dict['gridpdf_suffix'] = 'grid.pdf'
name_dict['Stackname'] = 'Stacking_Masked_MasterGrid.pdf'
name_dict['Stackname_nomask'] = 'Stacking_MasterGrid.pdf'
name_dict['Average_Bin_Value'] = 'Average_R23_O32_Values.tbl'
name_dict['temp_metallicity_pdf'] = '_Temp_Composite_Metallicity.pdf'

bian_coeff = [-0.32293, 7.2954, -54.8284, 138.0430]


def read_fitsfiles(fits_file_path):
    fits_data, header = fits.getdata(fits_file_path, header=True)
    wave = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
    dispersion = header['CDELT1']

    # Initialize Dictionary
    fits_dict = dict()
    # Fill in Dictionary
    fits_dict['fits_data'] = fits_data
    fits_dict['header'] = header
    fits_dict['wave'] = wave
    fits_dict['dispersion'] = dispersion

    return fits_dict


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
    log.info(f"exclude flag: {np.where(exclude_flag == 1)[0]}")

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

    indv_prop_tab = Table(table_dict)

    # We can create two different kinds of tables here of the R23_032 data (det3)
    # used to be get_det3_table.tbl
    indv_prop_file = join(fitspath, filename_dict['indv_prop'])
    if not exists(indv_prop_file):
        log.info(f"Writing: {indv_prop_file}")
    else:
        log.info(f"Overwriting: {indv_prop_file}")
    asc.write(indv_prop_tab, indv_prop_file, format='fixed_width_two_line',
              overwrite=True)

    log.info("finished.")

    return individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, det3, data3
