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

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'


outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
voronoi = np.loadtxt(outfilevoronoi)
voronoi = voronoi.transpose()
xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]


RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
#field0 = range(len(rframe_files))



def Master_Stacking(wave, image2D, name, header, mask= None):
    pdf_pages = PdfPages(fitspath+name) #open pdf document 
    

    #image2D = np.array(image2DM)
    image2DM = np.nan_to_num(image2D) 
    
    if mask !=None:
        image2DM = np.ma.masked_array(image2DM, mask)       

    n_bins = np.int(np.max(voronoi[2]))
    
    outfile = name.replace('.pdf', '.fits')  #fits file name and initialization 
    if not exists(outfile):
        stack_2d = np.zeros((n_bins+1, len(wave)), dtype = np.float64)
    else:
        #print 'reading ', outfile
        stack_2d = fits.getdata(outfile)
        
    for rr in range(n_bins+1):
        #print rr, stack_2d.shape
        index= np.where(voronoi[2]== rr)[0]   
        subgrid= image2DM[index]

        if exists(outfile):
            Spect1D = stack_2d[rr]
        else:
            if mask != None:                         
                Spect1D = np.ma.mean(subgrid, axis=0)
                #Spect1D = np.nanmedian(subgrid, axis=0)
            else:
                Spect1D = np.nanmean(subgrid, axis=0)

        stack_2d[rr] = Spect1D
        #print 'Spect1D:' , Spect1D
        #print 'index : ', index

        #Compute number of spectra at a given wavelength
        a = ma.count(subgrid, axis=0)
                
        fig, ax = plt.subplots()
            
        #GridSpec
        gs  = GridSpec(6,1)
        ax1 = plt.subplot(gs[0,0])
        ax2  = plt.subplot(gs[1:,0])

        #ax1 Plot
        ax1.plot(wave, a, linewidth=0.5)
            
        
        #ax2 plot
        ax2.plot(wave, Spect1D, linewidth=0.5)
        ax2.set_xlabel('Wavelength')
        ax2.set_ylabel('Spectra 1D')
        ax2.minorticks_on()
        #ax.set_xlim(np.log10(xlim))
        ax2.set_ylim([0,np.nanmax(Spect1D)])

        ax2.legend(loc='upper left', numpoints=3)    
        
        #str0 = 'R23=%.1f O32=%.1f N=%i' % (R23_grid[rr], O32_grid[oo], len(index))
        #ax2.annotate(str0, (0.05,0.95), xycoords='axes fraction', ha='left', va='top', weight='bold')
        
        for x in xcoor: ax2.axvline(x=x, linewidth= 0.3, color= 'k')
        
        fig.set_size_inches(8,11)
        
        
        x1, x2  = 4200, 4500
        ind= np.where((wave>= x1) & (wave<= x2))[0]
        sig, med = np.nanstd(Spect1D[ind]), np.nanmedian(Spect1D[ind])
        y1, y2 = med-2.5*sig, np.nanmax(Spect1D[ind])
        axins = zoomed_inset_axes(ax2, 2, loc=6, bbox_to_anchor=[375, 425])  #loc=9 
        axins.plot(wave, Spect1D, linewidth=0.5)
        axins.set_xlim([x1, x2])
        axins.set_ylim([y1, y2])
        axins.minorticks_on()
                
        for x in xcoor: axins.axvline(x=x, linewidth= 0.3, color= 'r')
        
        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="k", ls='dashed', lw=0.5)

        fig.tight_layout()
        plt.draw()
        pdf_pages.savefig(fig)
    #endfor
    pdf_pages.close()

    #Writing fits file
    if not exists(outfile):
        print 'writing ', outfile
        header['CTYPE1'] = 'LINEAR'
        header['CTYPE2'] = 'LINEAR'
        header['CRPIX1'] =  1.00000
        header['CDELT2'] =  1.00000
        header['CRPIX2'] =  1.00000

        #i_x, i_y = np.where(stack_2d == 0)
        #print len(i_x), len(i_y)
        #stack_2d[i_x,i_y] = np.nan
        fits.writeto(outfile, stack_2d[0:n_bins], header, overwrite= True)
    



#print('Starting MasterGrid')
def run_Stacking_Master_mask():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    name = 'Stacking_Voronoi_masked.pdf'
    maskM = fits.getdata(fitspath+'/Results_Graphs_Arrays_Etc/Arrays/MastermaskArray.fits') 
    MasterStack0 = Master_Stacking(wavemaster, image2DM, name, header, mask=maskM)
    
def run_Stacking_Master():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    name = 'Stacking_Voronoi.pdf'
    MasterStack0 = Master_Stacking(wavemaster, image2DM, name, header)
    

#run_Stacking_Master_mask()
#run_Stacking_Master()

