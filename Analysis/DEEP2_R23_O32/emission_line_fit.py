
###THIS CODE IS THE FIRST ITERATION OF THE ZOOM_AND_GAUSS_GENERAL CODE
###NOT CALLED IN GENERAL FUNCTION


import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
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

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
voronoi = np.loadtxt(outfilevoronoi)
voronoi = voronoi.transpose()

stacking_vor= r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Stacking_Voronoi_masked_output.fits'

tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi.tbl'
asc_tab = asc.read(tab)

stack2D, header = fits.getdata(stacking_vor, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32,5006.84]

def gauss(x,xbar,s,a,c):
   
    return  a*np.exp(-(x-xbar)**2/s**2) + c 


x_idx = np.where((wave>=4800) & (wave<=5100))[0]   #wave[4300:4400    #np.arange(0,100,1)
x_idx2 = np.where((wave>=4950) & (wave<=5050))[0]
xbar = 5007. #emission line 
s=1.0
a= 3.0
c = 2.0

#y  = gauss(wave[x_idx], xbar,s, a, c)


def get_func():
 
    #f = gauss(x,s,a,c)

    y0 = stack2D[2]
    x0 = wave[x_idx2]

    med0 = np.median(y0[x_idx2])
    max0 = np.max(y0[x_idx]) - med0
    print 'med0:', med0
    print 'max0:', max0
    p0 = [5006.8, 1.0, max0, med0] #must have some reasonable values
    o1, o2 = curve_fit(gauss, x0, y0[x_idx2], p0=p0,sigma=None) #verbose= True)
    #plt.plot(x,f)
    print o1
    print o2
    plt.plot(wave, y0,'r')
    plt.axhline(y=med0, linewidth= 0.5, color= 'k')
    plt.plot(x0,gauss(x0,*o1),'b--')
    plt.xlim(4800,5100)
    plt.show()



def all_fit():
    for rr in xcoor:
        name = 'Stacking_Voronoi_Zoomed_new.pdf'
        y0 = stack2D[2]
        x0 = wave[x_idx2]
        pdf_pages = PdfPages(fitspath+name)
        
        med0 = np.median(y0[x_idx2])
        max0 = np.max(y0[x_idx]) - med0
        print 'med0:', med0
        print 'max0:', max0
        p0 = [rr, 1.0, max0, med0] #must have some reasonable values
        o1, o2 = curve_fit(gauss, x0, y0[x_idx2], p0=p0,sigma=None) #verbose= True)
        #plt.plot(x,f)
        print o1
        print o2
        plt.plot(wave, y0,'r')
        plt.axhline(y=med0, linewidth= 0.5, color= 'k')
        plt.plot(x0,gauss(x0,*o1),'b--')
        x_min= rr-100
        x_max= rr+100
        plt.xlim(4800,5100)
        plt.plot.set_size_inches(11,8)
        plt.draw()
        plt.savefig(pdf_pages, format='pdf')
        
    

                                             
   
    pdf_pages.close()
        
