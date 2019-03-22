
###ORIGINAL FILE THAT RMS FUNCTION WAS WRITTEN IN
###NOW THE RMS FUNCTION IS IN THE ZOOM_AND_GAUSS_GENERAL CODE
###THIS FUNCTION IS NOT RUN



import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ



fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/combined_flux_table.tbl'
asc_tab = asc.read(tab)


stacking_vor= r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Stacking_Voronoi_masked_output.fits'
stack2D, header = fits.getdata(stacking_vor, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
dispersion = header['CDELT1']
Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')

lambda0 =[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  #11

lineflag = np.zeros(len(wave))
for ii in lambda0:   
    idx = np.where(np.absolute(wave - ii)<=5)[0]
    if len(idx) > 0:
        lineflag[idx] = 1

#Noise = sigma(f) * spectral dispersion(wavelast-wavefirst)* sqrt(number of pixals along the distribution)

#Order: ['OII_3727','H_9','NeIII','HeI','HEPSIL', 'HDELTA','HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007'


def rms_func(lambda_in, line_name):
    x_idx = np.where((wave-lambda_in)<=100 & (lineflag==0))[0]
    ini_sig = np.zeros(Spect_1D.shape[0])
    for rr in range(Spect_1D.shape[0]):
        y0 = stack2D[rr]
        sigma = np.std(y0[x_idx])

       
        pix =  5* sigma_array[rr]* dispersion 
        s_pix = np.sqrt(pix)

        ini_sig[rr]= s_pix * sigma * dispersion 
    return ini_sig
       

    


