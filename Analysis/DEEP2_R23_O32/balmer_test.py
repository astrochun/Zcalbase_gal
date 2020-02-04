##Testing Balmer Emission Line Fits###


import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc

import zoom_and_gauss_general


fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/balmer_test/'

'''lambda0 = [3797.90, 3835.38, 3889.05, 3970.07, 4101.73, 4340.46, 4861.32]  #, '4101.73, 4363.21, 4861.32']
        #[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Balmer', 'Balmer', 'Balmer', 'Balmer', 'Balmer', 'Balmer', 'Balmer']

line_name = ['H_10', 'H_9','Hzeta', 'HEPSIL', 'HDELTA', 'Hgamma', 'HBETA']'''


lambda0 = [4861.32]  #, '4101.73, 4363.21, 4861.32']
        #[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Balmer']

line_name = ['HBETA']

y_correction=''

#Double
dataset = 'Double_Bin'
Stack_name = 'Stacking'+dataset+'_output.pdf'
Stack_name = Stack_name.replace('.pdf', '.fits')
outfile_grid = fitspath + Stack_name
print outfile_grid
stack2D, header = fits.getdata(outfile_grid, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
#Spect_1D = fits.getdata(outfile_grid)
dispersion = header['CDELT1']
binning_avg_asc = fitspath+'/'+dataset+'_binning_averages.tbl'

lineflag = np.zeros(len(wave))
for ii in lambda0:   
    idx = np.where(np.absolute(wave - ii)<=5)[0]
    if len(idx) > 0:
        lineflag[idx] = 1
#Option to change: Constants used as initial guesses for gaussian fit
s=1.0
a= 1.0
c = 2.0
        
s1= 1.3
a1= 4.7
s2 = 10.0
a2 = -2.0

zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name, y_correction, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)



'''
#Voronoi14
dataset = 'Voronoi14'
Stack_name = 'Stacking'+dataset+'output.pdf'
Stack_name = Stack_name.replace('.pdf', '.fits')
outfile_grid = fitspath + Stack_name
print outfile_grid
stack2D, header = fits.getdata(outfile_grid, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
#Spect_1D = fits.getdata(outfile_grid)
dispersion = header['CDELT1']
binning_avg_asc = fitspath+'/'+dataset+'binning_averages.tbl'

lineflag = np.zeros(len(wave))
for ii in lambda0:   
    idx = np.where(np.absolute(wave - ii)<=5)[0]
    if len(idx) > 0:
        lineflag[idx] = 1
#Option to change: Constants used as initial guesses for gaussian fit
s=1.0
a= 1.0
c = 1
        
s1=-0.3
a1= 4.7
s2 = 1
a2 = -1.8

zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)






#Voronoi20
dataset = 'Voronoi20'
Stack_name = 'Stacking'+dataset+'output.pdf'
Stack_name = Stack_name.replace('.pdf', '.fits')
outfile_grid = fitspath + Stack_name
print outfile_grid
stack2D, header = fits.getdata(outfile_grid, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
#Spect_1D = fits.getdata(outfile_grid)
dispersion = header['CDELT1']
binning_avg_asc = fitspath+'/'+dataset+'binning_averages.tbl'

lineflag = np.zeros(len(wave))
for ii in lambda0:   
    idx = np.where(np.absolute(wave - ii)<=5)[0]
    if len(idx) > 0:
        lineflag[idx] = 1
#Option to change: Constants used as initial guesses for gaussian fit
s=1.0
a= 1.0
c = 1
        
s1=-0.3
a1= 4.7
s2 = 1
a2 = -1.8

zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)
'''



###Plotting Balmer###


