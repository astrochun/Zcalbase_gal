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

#If 'x0' is infeasible error occurs, check the para_bound values to make sure the expected values are within the range set up upper and lower limits. 

RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
#eventually I need to incoorparate these values into the correct function but right now they are global variables 
s=1.0
a= 1.0
c = 2.0
s1= 1.3
a1= 1.5
s2 = 1
a2 = 1.8

x_signal = [4363.21,4958.91,5006.84, 3835.38, 3970.07, 4101.73, 4340.46, 4861.32]


def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth

def gauss(x,xbar,s,a,c):
   
    return  a*np.exp(-(x-xbar)**2/s**2) + c 

def double_gauss(x, xbar, s1, a1, c, s2, a2): 

    return a1*np.exp(-(x-xbar)**2/s1**2) + c + a2*np.exp(-(x-xbar)**2/s2**2)

def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    con1 = 3728.91/3726.16
    return a1*np.exp(-(x-xbar)**2/s1**2) + c + a2*np.exp(-(x-(xbar*con1))**2/s2**2) 

def get_gaussian_fit(working_wave,x0, y0, y_norm, x_idx,keyword):

    med0 = np.median(y_norm[x_idx])
    max0 = np.max(y_norm[x_idx]) - med0
    print med0, max0 
    if keyword == 'Single': 
        p0 = [working_wave, 1.0, max0, med0] #must have some reasonable values
        print p0  #p0 = [working_wave, sigma, amplitude (max0), constant (med0)]

        para_bound = ((working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0)),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0)))
        o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)
        print o1

    if keyword == 'Balmer':
        p0 = [working_wave, 1.0, max0, med0,s2, -0.1*max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0,0)
        print para_bound
        print working_wave,  min(x0[x_idx]), max(x0[x_idx]), 

        o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)
    

    if keyword == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75*max0, med0, 1.0, max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, 0.0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0, 100.0)
       

        o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)

    return o1, med0, max0
  

#Plotting Zoomed 
def zoom_gauss_plot(working_wave,keyword = '',outpdf=''):
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
    if outpdf == '':
        name = 'Stacking_Voronoi_Zoomed_Gauss_generalexperiment.pdf'
        outpdf = fitspath + name

    pdf_pages = PdfPages(outpdf)
    nrows = 4
    ncols = 4
    x_idx = np.where((wave>=(working_wave-100)) & (wave<=(working_wave+100)))[0] 
    x0 = wave#[x_idx]
    scalefact = 1e-17
    
    for rr in range(Spect_1D.shape[0]):
        y0 = stack2D[rr]
        y_norm = y0/scalefact

        row = rr / nrows % ncols
        col = rr % ncols
        print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1 = working_wave-100
        x2  = working_wave+100

        y_smooth = movingaverage_box1D(Spect_1D[rr]/scalefact, 5, boundary='extend')
      

        o1, med0, max0  = get_gaussian_fit(working_wave,x0, y0, y_norm, x_idx, keyword)
     

        #Calculating Flux: Signal Line Fit
        x_sigsnip = np.where((np.abs((x0-working_wave))/np.abs(o1[1]))<=2.5)[0]
        if keyword == 'Single':
            gauss0=gauss(x0[x_sigsnip],*o1)
        if keyword == 'Balmer':
            gauss0 = double_gauss(x0[x_sigsnip], *o1)
        if keyword == 'Oxy2':
            gauss0 = oxy2_gauss(x0[x_sigsnip], *o1)
        dx = x0[2]-x0[1]  #0.16866575
        flux_g = np.sum((gauss0-o1[3])*dx)    #flux from gaussian distribution 
        
        flux_s= np.sum((y_norm[x_sigsnip]-o1[3])*dx)    #flux from snipping method (spectral flux)where snip off sigma >2.5
        if rr == 0: print 'o1', o1, flux_g, flux_s, x_sigsnip


        #Plotting
        emis= t_ax.plot(wave, y_smooth,'k', linewidth=0.75, label= 'Emission')
        if keyword == 'Oxy2':
            t_ax.set_xlim([x1+45,x2-45])
        else:
            t_ax.set_xlim([x1,x2])

        if keyword == 'Single':
            gauss0=gauss(x0,*o1)
        if keyword == 'Balmer':
            gauss0= double_gauss(x0,*o1)
        if keyword == 'Oxy2':
            gauss0 = oxy2_gauss(x0,*o1)


        t_ax.plot(x0,gauss0, 'b--', linewidth= 0.65, label= 'Gauss Fit')
        t_ax.legend(loc=4, ncol=2, fontsize = 6)
        
        txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 += 'R_23%.3f O32 %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])  #$\overline{x}$:$\overline{y}$:
        txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr]) 
        txt0 += 'Median: %.3f Sigma: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
        txt0 += 'Flux_G: %.3f Flux_S %.3f' %(flux_g, flux_s)
        
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')
        for x in x_signal: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        #if row == nrows-1: t_ax.set_xlabel('Wavelength')

        if col == 0:
             t_ax.set_ylabel('Spect_1D')
        t_ax.set_yticklabels([])  #sets y-tick labels 
    
        if rr == Spect_1D.shape[0]-1 and rr % (nrows*ncols) != nrows*ncols-1:
            for jj in range(col+1, ncols): ax_arr[row,jj].axis('off')
            for kk in range(row+1, nrows):
                for zz in range(ncols): ax_arr[kk,zz].axis('off')

            xlabels = t_ax.get_xticklabels()
            print 'xlabels : ', xlabels
            for yy in range(rr,rr-ncols,-1):
                y_col = yy % ncols
                y_row = yy / nrows % ncols
                ax_arr[y_row,y_col].set_xlabel('Wavelength')
                ax_arr[y_row,y_col].set_xticklabels(xlabels)
                
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97, hspace=0.05)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(8,8)
            plt.draw()
            fig.savefig(pdf_pages, format='pdf')
        
    

                                             
    #endfor
    pdf_pages.close()
    print 'Done!'



def zm_general():
    #x_signal = [4363.21,4958.91,5006.84]
    #x_balmer = [3835.38, 3970.07, 4101.73, 4340.46, 4861.32] #3797.90
    s=1.0
    a= 1.0

    c = 1
    s1=-0.3
    a1= 4.7
    s2 = 1
    a2 = -1.8
    key = 'Oxy2' #'Balmer' # 'Single'
    '''for aa in x_signal:
        outpdf = fitspath+'Stacking_Voronoi_Zoomed_Gauss_'+str(np.int(aa))+'.pdf'
        print outpdf
        m0 = zoom_gauss_plot(aa, keyword= key, outpdf=outpdf)'''

    '''for aa in x_balmer:
        outpdf = fitspath+'Stacking_Voronoi_Zoomed_Balmer_'+str(np.int(aa))+'.pdf'
        print outpdf
        m0 = zoom_gauss_plot(aa, keyword=key, outpdf=outpdf)'''

    if key == 'Oxy2':
        working_wave = 3726.16 #3727.5
        outpdf = fitspath+'Stacking_Voronoi_Zoomed_OxygenII_'+str(np.int(working_wave))+'.pdf'
        print outpdf
        m0 = zoom_gauss_plot(working_wave, keyword=key, outpdf=outpdf)


#erg/s/sm^2/A
#plot residuals 
#learn how to integrate


