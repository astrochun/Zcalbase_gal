##In this file, I define the stacking code in a function that runs over the master grid
##Creates a pdf and a fits file
##fits file for his code is given the name of the PDF file 
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

#fitspath='/astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/'
fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'


outfile01 = 'Arrays_R23O32bin01MasterGrid.npz'
#outfile025 = '/astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/Arrays_R23O32bin025MasterGrid.npz'
outfile025 = fitspath + '/Results_Graphs_Arrays_Etc/Arrays_R23O32bin025MasterGrid.npz'
grid_data = np.load(outfile025)
xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]


RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits' 


def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth

def Master_Stacking(wave, image2D, name, header, mask= None):
    pdf_pages = PdfPages(fitspath+name) #open pdf document 
    
    R23_grid = grid_data['R23_grid']
    O32_grid = grid_data['O32_grid']
    #print len(grid_data['T_arr'][:,0])
    #print len(grid_data['T_arr'][0,:])
    #print len(R23_grid), len(O32_grid)

    #image2D = np.array(image2DM)
    image2DM = np.nan_to_num(image2D) 
    
    if mask !=None:
        image2DM = np.ma.masked_array(image2DM, mask)       

    outfile = name.replace('.pdf', '.fits')  #fits file name and initialization 
    if not exists(outfile):
        stack_2d = np.zeros((len(R23_grid)*len(O32_grid), len(wave)), dtype=np.float64)
    else:
        print 'reading ', outfile
        stack_2d = fits.getdata(outfile)
        
    count = 0
    for rr in range(len(R23_grid)):
        for oo in range(len(O32_grid)):
            index= grid_data['T_arr'][rr,oo]
            if len(index) >10:   
                subgrid= image2DM[index]

                if exists(outfile):
                    Spect1D = stack_2d[count]
                else:
                    if mask != None:                         
                        Spect1D = np.ma.mean(subgrid, axis=0)
                        #Spect1D = np.nanmedian(subgrid, axis=0)
                    else:
                        Spect1D = np.nanmean(subgrid, axis=0)

                    stack_2d[count] = Spect1D
                print "stack_2d",[rr], Spect1D
                if rr ==0 and oo == 0: print index

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

                str0 = 'R23=%.1f O32=%.1f N=%i' % (R23_grid[rr], O32_grid[oo], len(index))
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

                # draw a bbox of the region of the inset axes in the parent axes and
                # connecting lines between the bbox and the inset axes area
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
    #if not exists(outfile):
    print 'writing ', outfile
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] =  1.00000
    header['CDELT2'] =  1.00000
    header['CRPIX2'] =  1.00000

        #i_x, i_y = np.where(stack_2d == 0)
        #print len(i_x), len(i_y)
        #stack_2d[i_x,i_y] = np.nan
    fits.writeto(outfile, stack_2d[0:count], header, overwrite= True)
    

print('Starting MasterGrid')
def run_Stacking_Master_mask():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    name = 'Stacking_Wave_vs_Spect1D_Masked_MasterGrid_Average_bin025.pdf'
    maskM = fits.getdata(fitspath+'/Results_Graphs_Arrays_Etc/Arrays/MastermaskArray.fits') 
    MasterStack0 = Master_Stacking(wavemaster, image2DM, name, header, mask=maskM)
    
def run_Stacking_Master():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    name = 'Stacking_Wave_vs_Spect1D_MasterGrid_Average_bin025.pdf'
    MasterStack0 = Master_Stacking(wavemaster, image2DM, name, header)
    

#run_Stacking_Master_mask()
#run_Stacking_Master()

#Plotting Zoomed in on 4363
def zoom_plot_4363():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'/Results_Graphs_Arrays_Etc/Arrays/Stacking_Wave_vs_Spect1D_Masked_MasterGrid_Average_bin025.fits')
    name = 'Stacking_Zoomed_4363.pdf'
    pdf_pages = PdfPages(fitspath+name)
    nrows = 4
    ncols = 4
    R23_grid = grid_data['R23_grid']
    O32_grid = grid_data['O32_grid']
    #index= grid_data['T_arr'][rr,rr]
    print Spect_1D.shape
    
    for rr in range(Spect_1D.shape[0]):

        row = rr / nrows % ncols
        col = rr % ncols
        print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1, x2  = 4200, 4500

        y_smooth = movingaverage_box1D(Spect_1D[rr], 5, boundary='extend')

        t_ax.plot(wave, y_smooth, linewidth=0.5)
        t_ax.set_xlim([x1,x2])
        #t_ax.set_xlabel('Wavelength')
        #t_ax.set_ylabel('Spectra 1D')
        for x in xcoor: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

        #txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 = r'R_23%.3f O32 %.3f\n' % (grid_data['R23_grid'][rr], grid_data['O32_grid'][rr])  #$\overline{x}$:$\overline{y}$:
        #txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr])
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')

        str0 = 'R23=%.1f O32=%.1f ' % (R23_grid[rr], O32_grid[rr])# N=%i  len(index))
        t_ax.annotate(str0, (0.05,0.95), xycoords='axes fraction', ha='left', va='top', weight='bold')
               

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        if col == 0:
             t_ax.set_ylabel('Spect_1D')
        else: t_ax.set_yticklabels([])
    
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

