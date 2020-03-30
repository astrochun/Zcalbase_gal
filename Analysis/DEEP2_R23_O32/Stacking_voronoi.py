
###THIS CODE DOES THE STACKING FOR ALL VORONOI AND (CURRENTLY) THE DOUBLE BINNING METHOD
###CURRENTLY RE-WRITING TO TEST THAT VORONOI AND GRID STACKING ARE THE SAME


import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
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
from astropy.table import vstack, hstack


xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32,5006.84]

from . import general 

RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'


def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth


def Master_Stacking(fitspath, dataset, voronoi_data, asc_table1, wave, image2D, name, header, mask= None):
    pdf_pages = PdfPages(fitspath+name) #open pdf document 
    
    voronoi_data = asc.read(voronoi_data)
    asc_tab = asc.read(asc_table1)
    print('vornoni_data_file:', voronoi_data)
    print('asc_tab:,', asc_tab)
    
    image2DM = np.nan_to_num(image2D) 
    
    if mask !=None:
        image2DM = np.ma.masked_array(image2DM, mask)       

    last_column = voronoi_data['N_bin'].data
    n_bins = np.max(last_column)    
    outfile = name.replace('.pdf', '.fits')  #fits file name and initialization 
    if not exists(outfile):
        print("File not found : "+outfile+".  Initiating!!!")
        stack_2d = np.zeros((n_bins+1, len(wave)), dtype = np.float64) 
    else:
        print("File found : "+outfile+".  Reading in!!!")
        stack_2d = fits.getdata(outfile)
        
    for rr in xrange(n_bins+1):  
        print(rr, stack_2d.shape)
        index= np.where(last_column== rr)[0]   
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
        y_smooth = movingaverage_box1D(Spect1D, 5, boundary='extend')
        ax2.plot(wave, y_smooth, linewidth=0.5)
        ax2.set_xlabel('Wavelength')
        ax2.set_ylabel('Spectra 1D')
        ax2.minorticks_on()
        ax2.set_ylim([0,np.nanmax(Spect1D)])

        ax2.legend(loc='upper left', numpoints=3)    
        
        for x in xcoor: ax2.axvline(x=x, linewidth= 0.3, color= 'k')

        if dataset == 'Double_Bin':
            txt0 = r'R23_start= %.3f  O32_median= %.3f  ' % (asc_tab['R23_value'][rr], asc_tab['O32_value'][rr]) + '\n'
            txt0 += 'ID: %.3f  N: %.3f\n' % (asc_tab['ID'][rr], asc_tab['area'][rr]) 
        else: 
            txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
            txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])
            txt0 += 'Scale: %.3f N: %.3f\n' % (asc_tab['scale'][rr], asc_tab['area'][rr]) 
        
       
        ax2.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '12')
       
        
        fig.set_size_inches(8,11)
        fig.tight_layout()
        plt.draw()
        pdf_pages.savefig(fig)
    #endfor
    pdf_pages.close()

    #Writing fits file
    #if not exists(outfile):
    print('writing ', outfile)
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] =  1.00000
    header['CDELT2'] =  1.00000
    header['CRPIX2'] =  1.00000

    T_fits_name = fitspath + outfile
    fits.writeto(T_fits_name, stack_2d, header, overwrite= True)

   
    fig.clear()
#Function that runs Master_Stacking and calls necessary inputs (including mask)
def run_Stacking_Master_mask(fitspath_ini, dataset, fitspath,voronoi_data, det3,asc_table1, Stack_name):
    
   
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    image2D = image2DM[det3]
    print(len(image2D))
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    maskM = fits.getdata(fitspath_ini+'/Results_Graphs_Arrays_Etc/Arrays/MastermaskArray.fits')
    mask= maskM[det3]
    MasterStack0 = Master_Stacking(fitspath, dataset, voronoi_data, asc_table1, wavemaster, image2D, Stack_name, header, mask=mask)



#Function that runs Master_Stacking and calls necessary inputs (without mask)    
def run_Stacking_Master(fitspath,dataset, voronoi_data, Stack_name):
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    MasterStack0 = Master_Stacking(fitspath, dataset, voronoi_data, asc_table1, wavemaster, image2DM, Stack_name, header, mask=None)
    



