# Testing Balmer Emission Line Fits

import numpy as np
from astropy.io import fits

import zoom_and_gauss_general

fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/balmer_test/'

lambda0 = [4861.32]  #, '4101.73, 4363.21, 4861.32']

line_type = ['Balmer']

line_name = ['HBETA']

y_correction = ''

# Option to change: Constants used as initial guesses for gaussian fit
s = 1.0
a = 1.0
c = 2.0

s1 = 1.3
a1 = 4.7
s2 = 10.0
a2 = -2.0


def balmer_fitting_test(dataset):
    # Double
    if dataset == 'Double_Bin':
        Stack_name = 'Stacking'+dataset+'_output.pdf'
        Stack_name = Stack_name.replace('.pdf', '.fits')
        outfile_grid = fitspath + Stack_name
        print outfile_grid
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        # Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']
        binning_avg_asc = fitspath+'/'+dataset+'_binning_averages.tbl'

        lineflag = np.zeros(len(wave))
        for ii in lambda0:
            idx = np.where(np.absolute(wave - ii) <= 5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1

    # Voronoi14
    if dataset == 'Voronoi14':
        Stack_name = 'Stacking' + dataset + 'output.pdf'
        Stack_name = Stack_name.replace('.pdf', '.fits')
        outfile_grid = fitspath + Stack_name
        print outfile_grid
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
        # Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']
        binning_avg_asc = fitspath + '/' + dataset + 'binning_averages.tbl'

        lineflag = np.zeros(len(wave))
        for ii in lambda0:
            idx = np.where(np.absolute(wave - ii) <= 5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1

    # Voronoi20
    if dataset == 'Voronoi20':
        Stack_name = 'Stacking' + dataset + 'output.pdf'
        Stack_name = Stack_name.replace('.pdf', '.fits')
        outfile_grid = fitspath + Stack_name
        print outfile_grid
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
        # Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']
        binning_avg_asc = fitspath + '/' + dataset + 'binning_averages.tbl'

        lineflag = np.zeros(len(wave))
        for ii in lambda0:
            idx = np.where(np.absolute(wave - ii) <= 5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1

    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave,
                                      lineflag, dispersion, lambda0, line_type, line_name,
                                      y_correction, s, a, c, s1, a1, s2, a2, tab=binning_avg_asc)
