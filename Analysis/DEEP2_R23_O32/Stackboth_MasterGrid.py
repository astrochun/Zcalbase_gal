
###THIS FUNCTION DOES THE STACKING FOR THE ORIGINAL AND SIGNAL BINNING METHODS
###CURRENTLY RE-WRITING TO TEST THAT VORONOI AND GRID STACKING ARE THE SAME

##In this file, I define the stacking code in a function that runs over the master grid
##Creates a pdf and a fits file
##fits file for this code is given the name of the PDF file 
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

import general

xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]


RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits' 


def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth

def Master_Stacking(fitspath,dataset, wave, grid_data_file, image2D, name, header, mask= None):
    pdf_pages = PdfPages(fitspath+name) #open pdf document 

    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3= general.get_det3()


    grid_data = np.load(grid_data_file)
    R23_grid = grid_data['R23_grid']
    O32_grid = grid_data['O32_grid']
    
    image2DM = np.nan_to_num(image2D[det3]) 
    
    print 'Mask[det3]', mask[det3]
    if type(mask) != type(None):
        image2DM = np.ma.masked_array(image2DM, mask[det3])
        
    outfile = name.replace('.pdf', '.fits')  #fits file name and initialization 
    if not exists(outfile):
        stack_2d = np.zeros((len(R23_grid)*len(O32_grid), len(wave)), dtype=np.float64)
    else:
        #print 'reading ', outfile
        stack_2d = fits.getdata(outfile)

    avg_R23 = np.zeros(len(R23_grid)*len(O32_grid))
    avg_O32 = np.zeros(len(R23_grid)*len(O32_grid))
    R23_node = np.zeros(len(R23_grid)*len(O32_grid))
    O32_node = np.zeros(len(R23_grid)*len(O32_grid))
    N_gal = np.zeros(len(R23_grid)*len(O32_grid))

    n_N = R23_grid.shape[0]
    if dataset == 'n_Bins' or dataset == 'Double_Bin': n_M = R23_grid.shape[1]
    else: n_M = O32_grid.shape[0]
    count = 0
    
    for rr in range(n_N):
        for oo in range(n_M):
            print rr,oo
            index= grid_data['T_arr'][rr,oo]
            #print rr, oo, 'index:', index, index.dtype
            if len(index) >10:
                R23_node[count] = R23_grid[rr,oo] 
                O32_node[count] = O32_grid[rr,oo]
                avg_R23[count]  = np.average(np.log10(R23)[index])
                avg_O32[count]  = np.average(np.log10(O32)[index])
                N_gal[count] = len(index)
                subgrid= image2DM[index]


                print('R23:', R23_node[count],'O32:', O32_node[count], 'avg_R23', avg_R23[count], 'avg_O32', avg_O32[count])

                if exists(outfile):
                    Spect1D = stack_2d[count]
                else:
                    if type(mask) != type(None):                         
                        Spect1D = np.ma.mean(subgrid, axis=0)
                
                    else:
                        Spect1D = np.nanmean(subgrid, axis=0)

                    stack_2d[count] = Spect1D

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
                ax2.set_ylim([0,np.nanmax(Spect1D)])

                ax2.legend(loc='upper left', numpoints=3)    

                str0 = 'R23=%.1f O32=%.1f N=%i' % (R23_grid[rr,oo], O32_grid[rr,oo], len(index))
                ax2.annotate(str0, (0.05,0.95), xycoords='axes fraction', ha='left', va='top', weight='bold')
               
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

                mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="k", ls='dashed', lw=0.5)

                fig.tight_layout()
                plt.draw()
                pdf_pages.savefig(fig)


                count +=1
            #endif else: count +=1
        #endfor
    #endfor
    pdf_pages.close()

    #Writing fits file
    print 'writing ', outfile
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] =  1.00000
    header['CDELT2'] =  1.00000
    header['CRPIX2'] =  1.00000

    fits.writeto(fitspath+outfile, stack_2d[0:count], header, overwrite= True)

    #Writing Ascii Tables and Fits Tables
    out_ascii = fitspath+'/'+dataset+'binning_averages.tbl'

    ID = np.arange(0,len(R23_node), 1, dtype = int)
    n=  ('ID','xnode', 'ynode', 'xBar', 'yBar', 'area')
    tab0 = Table([ID, R23_node, O32_node,avg_R23,avg_O32, N_gal], names=n)
    asc.write(tab0[0:count], out_ascii, format='fixed_width_two_line')

    fig.clear()    
def run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, name,grid_data_file):
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    maskM = fits.getdata(fitspath_ini+'/Results_Graphs_Arrays_Etc/Arrays/MastermaskArray.fits') 
    MasterStack0 = Master_Stacking(fitspath,dataset, wavemaster, grid_data_file,image2DM, name, header, mask=maskM)
    
def run_Stacking_Master(fitspath, name,grid_data_file):
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    MasterStack0 = Master_Stacking(fitspath, wavemaster,grid_data_file, image2DM, name, header)
    
