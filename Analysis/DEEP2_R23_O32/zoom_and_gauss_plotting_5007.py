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

RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'



def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth

def gauss(x,xbar,s,a,c):
   
    return  a*np.exp(-(x-xbar)**2/s**2) + c 


x_idx = np.where((wave>=4800) & (wave<=5100))[0]   #wave[4300:4400    #np.arange(0,100,1)
x_idx2 = np.where((wave>=4950) & (wave<=5050))[0]
xbar = 5007. #emission line 
s=1.0
a= 3.0
c = 2.0

#y  = gauss(wave[x_idx], xbar,s, a, c)

#y0 = stack2D[2]
x0 = wave[x_idx2]

def get_func(y0):
 
    #f = gauss(x,s,a,c)

   

    med0 = np.median(y0[x_idx2])
    max0 = np.max(y0[x_idx]) - med0
    #print 'med0:', med0
    #print 'max0:', max0
    p0 = [5006.8, 1.0, max0, med0] #must have some reasonable values
    #param_bounds = ((0, xval[np.argmax(yval)]-3, 0, -0.1*np.max(yval), xval[np.argmax(yval)]-3, 0, med0-0.05*np.abs(med0)),
        #(1e-15/1e-17, xval[np.argmax(yval)]+3, 10, 0, xval[np.argmax(yval)]+3, 10, med0+0.05*np.abs(med0)))
    #(x_min, )(x_max, )
    o1, o2 = curve_fit(gauss, x0, y0[x_idx2], p0=p0,sigma=None) #verbose= True)
    #plt.plot(x,f)
    #print o1
    #print o2
    #plt.plot(wave, y0,'r')
    #plt.axhline(y=med0, linewidth= 0.5, color= 'k')
    #plt.plot(x0,gauss(x0,*o1),'b--')
    #plt.xlim(4800,5100)
    #print o1
    return o1, med0, max0
    #plt.show()


#Plotting Zoomed in on 5007    
def zoom_gauss_plot_5007():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
    name = 'Stacking_Voronoi_Zoomed_Gauss_5007.pdf'
    pdf_pages = PdfPages(fitspath+name)
    nrows = 4
    ncols = 4
    
    for rr in range(Spect_1D.shape[0]):
        y0 = stack2D[rr]
        x0 = wave[x_idx2]
        row = rr / nrows % ncols
        col = rr % ncols
        print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1, x2  = 4950, 5050

        y_smooth = movingaverage_box1D(Spect_1D[rr], 5, boundary='extend')
        
        o1, med0, max0  = get_func(y0)

        #Calculating Flux
        x_sigsnip = np.where(((x0-5007)/o1[1])<=2.5)[0]
        gauss0=gauss(x0[x_sigsnip],*o1)
        dx = x0[2]-x0[1]  #0.16866575

        flux_g = np.sum((gauss0-med0)*dx)    #flux from gaussian distribution 
        ind = np.where((x0-5007)/o1[1] <=2.5)[0]
        flux_s= np.sum((y0-med0)[ind]*dx)    #flux from snipping method (spectral flux)where snip off sigma >2.5


        #Plotting
        emis= t_ax.plot(wave, y_smooth,'k', linewidth=0.75, label= 'Emission')
        t_ax.set_xlim([x1,x2])

        gauss0=gauss(x0,*o1)
        t_ax.plot(x0,gauss0, 'b--', linewidth= 0.65, label= 'Gauss Fit')
        t_ax.legend(loc=4, ncol=2, fontsize = 6)
        
        txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 += 'R_23%.3f O32 %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])  #$\overline{x}$:$\overline{y}$:
        txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr]) 
        txt0 += 'Median: %.3f Sigma: %.3f  Norm: %.3f'% (med0/1e-18, o1[1], max0/1e-17) + '\n'
        txt0 += 'Flux_G 5007: %.3f Flux_S 5007 %.3f' %(flux_g/1e-17, flux_s/1e-17)
        
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        if col == col:
             t_ax.set_ylabel('Spect_1D')
        #else: t_ax.set_yticklabels([])
    
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.2, right=0.98, bottom=0.06, top=0.95,
                            hspace=0.0)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(11,8)
            plt.draw()
            fig.savefig(pdf_pages, format='pdf')
        
    

                                             
    #endfor
    pdf_pages.close()
    print 'Done!'

#Calculating the Flux of each fit
def calcul_flux():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')

    for rr in range(Spect_1D.shape[0]):
        y0 = stack2D[rr]
        x0 = wave[x_idx2]
        o1, med0, max0  = get_func(y0)

        x1 = np.where(((x0-5007)/o1[1])<=2.5)[0]
        gauss0=gauss(x0[x1],*o1)
        dx = x0[2]-x0[1]  #0.16866575

        flux_g = np.sum((gauss0-med0)*dx)    #flux from gaussian distribution 
        ind = np.where((x0-5007)/o1[1] <=2.5)[0]
        flux_s= np.sum((y0-med0)[ind]*dx)    #flux from snipping method where snip off sigma >2.5
        print flux_g, flux_s











""" flux = sum(y0-med0)[index] dx
dx is the step size in the data
np.where ((x0-5007))/sigma <= 2.5) """
