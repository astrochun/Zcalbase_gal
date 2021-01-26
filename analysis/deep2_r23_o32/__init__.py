from astropy.io import fits
import numpy as np

# Dictionary of pdf and other files specific to the Zcalbase_gal project
name_dict = dict()

name_dict['gridnpz_suffix'] = 'grid.npz'
name_dict['gridpdf_suffix'] = 'grid.pdf'
name_dict['Stackname'] = 'Stacking_Masked_MasterGrid.pdf'
name_dict['Stackname_nomask'] = 'Stacking_MasterGrid.pdf'
name_dict['Average_Bin_Value'] = 'Average_R23_O32_Values.tbl'
name_dict['temp_metallicity_pdf'] = '_Temp_Composite_Metallicity.pdf'


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