
"""
THIS FUNCTION DOES THE STACKING FOR THE ALL BINNING METHODS EXCEPT VORNOI

In this file, I define the stacking code in a function that runs over the master grid
Creates a pdf and a fits file
fits file for this code is given the name of the PDF file

Emission lines in spectrum (not all being used currently in study) See MSC for subset
[3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46,
4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists, join
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from astropy.convolution import Box1DKernel, convolve

from . import general
from Metallicity_Stack_Commons import lambda0
from Metallicity_Stack_Commons.column_names import filename_dict

RestframeMaster = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'


def movingaverage_box1d(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth


def master_stacking(fitspath, dataset, wave, grid_data_file, image2D, name, header, mask= None):
    """
    Purpose
    Function stacks all spectra in a given bin and produces tables of properties of that bin

    Parameters
    fitspath -> save location of the current run
    dataset -> keyword used to define binning method
    wave -> spectrum of each galaxy
    grid_data_file -> npz file that holds the information from the binning process
    image2D -> spectra data
    name -> name of the outputted pdf file with graphs
    header -> header of the data file
    mask -> optional input used to mask the night sky lines if inputted (default: None)
    """
    pdf_pages = PdfPages(fitspath+name) # open pdf document

    individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3= general.get_det3(fitspath)


    grid_data = np.load(grid_data_file, allow_pickle=True) # This is the npz file
    R23_minimum = grid_data['R23_minimum']
    O32_minimum = grid_data['O32_minimum']

    image2DM = np.nan_to_num(image2D[det3])
    
    print('Mask[det3]', mask[det3])
    if type(mask) != type(None):
        image2DM = np.ma.masked_array(image2DM, mask[det3])
        
    outfile = join(fitspath, filename_dict['comp_spec'])
    if not exists(outfile):
        stack_2d = np.zeros((len(R23_minimum)*len(O32_minimum), len(wave)), dtype=np.float64)
    else:
        stack_2d = fits.getdata(outfile)
 
    avg_R23 = np.zeros(len(R23_minimum)*len(O32_minimum))    #Same as xBar
    avg_O32 = np.zeros(len(R23_minimum)*len(O32_minimum))    #Same as yBar
    R23_node = np.zeros(len(R23_minimum)*len(O32_minimum))    #Same as R23_minimum
    O32_node = np.zeros(len(R23_minimum)*len(O32_minimum))    #Same as O32_minimum
    R23_med = np.zeros(len(R23_minimum)*len(O32_minimum))    #median R23 value
    O32_med = np.zeros(len(R23_minimum)*len(O32_minimum))    #median O32 value
    N_gal = np.zeros(len(R23_minimum)*len(O32_minimum))     #Same as Number_inbin

    n_N = R23_minimum.shape[0]
    if dataset == 'n_Bins' or dataset == 'Double_Bin': n_M = R23_minimum.shape[1]
    else: n_M = O32_minimum.shape[0]
    count = 0
    
    for rr in range(n_N):
        for oo in range(n_M):
            print(rr,oo)
            index= grid_data['locator'][rr, oo]

            # calculating the average and minimum values of R23 and O32
            if len(index) >10:
                R23_node[count] = R23_minimum[rr, oo]
                O32_node[count] = O32_minimum[rr, oo]
                avg_R23[count] = np.average(R23[index])        #np.log10(R23)
                avg_O32[count] = np.average(O32[index]) #(np.log10(O32
                R23_med[count] = np.median(R23[index])   #np.log10(R23)
                O32_med[count] = np.median(O32[index])
                N_gal[count] = len(index)
                subgrid = image2DM[index]


                print('R23:', R23_node[count],'O32:', O32_node[count],
                      'avg_R23', avg_R23[count], 'avg_O32', avg_O32[count])

                if exists(outfile):
                    Spect1D = stack_2d[count]
                else:
                    if type(mask) != type(None):                         
                        Spect1D = np.ma.mean(subgrid, axis=0)
                
                    else:
                        Spect1D = np.nanmean(subgrid, axis=0)

                    stack_2d[count] = Spect1D

                # Compute number of spectra at a given wavelength
                a = ma.count(subgrid, axis=0)
                
                fig, ax = plt.subplots()

                # GridSpec
                gs = GridSpec(6,1)
                ax1 = plt.subplot(gs[0,0])
                ax2 = plt.subplot(gs[1:,0])

                # ax1 Plot
                ax1.plot(wave, a, linewidth=0.5)
                
                
                # ax2 plot
                ax2.plot(wave, Spect1D, linewidth=0.5)
                ax2.set_xlabel('Wavelength')
                ax2.set_ylabel('Spectra 1D')
                ax2.minorticks_on()
                ax2.set_ylim([0, np.nanmax(Spect1D)])

                ax2.legend(loc='upper left', numpoints=3)    

                str0 = 'R23=%.1f O32=%.1f N=%i' % (R23_minimum[rr, oo], O32_minimum[rr oo], len(index))
                ax2.annotate(str0, (0.05, 0.95), xycoords='axes fraction', ha='left', va='top', weight='bold')
               
                for x in lambda0:
                    ax2.axvline(x=x, linewidth= 0.3, color= 'k')

                fig.set_size_inches(8, 11)

                
                x1, x2  = 4200, 4500
                ind = np.where((wave >= x1) & (wave <= x2))[0]
                sig, med = np.nanstd(Spect1D[ind]), np.nanmedian(Spect1D[ind])
                y1, y2 = med-2.5*sig, np.nanmax(Spect1D[ind])
                axins = zoomed_inset_axes(ax2, 2, loc=6, bbox_to_anchor=[375, 425])  #loc=9 
                axins.plot(wave, Spect1D, linewidth=0.5)
                axins.set_xlim([x1, x2])
                axins.set_ylim([y1, y2])
                axins.minorticks_on()
                    
                for x in lambda0:
                    axins.axvline(x=x, linewidth= 0.3, color='r')

                mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="k", ls='dashed', lw=0.5)

                fig.tight_layout()
                plt.draw()
                pdf_pages.savefig(fig)


                count +=1
            # endif else: count +=1
        # endfor
    # endfor
    pdf_pages.close()

    # Writing fits file
    print('writing ', outfile)
    header['CTYPE1'] = 'LINEAR'
    header['CTYPE2'] = 'LINEAR'
    header['CRPIX1'] =  1.00000
    header['CDELT2'] =  1.00000
    header['CRPIX2'] =  1.00000

    fits.writeto(outfile, stack_2d[0:count], header, overwrite= True)

    # Writing Ascii Tables and Fits Tables
    out_ascii = fitspath+'/bin_info.tbl'   # used to be 'binning_averages.tbl'

    ID = np.arange(0,len(R23_node), 1, dtype = int)
    n=  ('bin_ID','logR23_min', 'logO32_min', 'logR23_avg', 'logO32_avg', 'logR23_med', 'logO32_med', 'N_stack')
    # for n_split xnode and ynode are the lowest values of the bin while xBar and yBar are the averages

    tab0 = Table([ID, R23_node, O32_node, avg_R23, avg_O32, R23_med, O32_med, N_gal], names=n)
    asc.write(tab0[0:count], out_ascii, format='fixed_width_two_line')

    fig.clear()


def run_stacking_master_mask(fitspath, fitspath_ini, dataset, name, grid_data_file):
    """
    Purpose
    Run function for above function if mask used
    """
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
    maskM = fits.getdata(fitspath_ini + '/Results_Graphs_Arrays_Etc/Arrays/MastermaskArray.fits')
    Master_Stacking(fitspath, dataset, wavemaster, grid_data_file, image2DM, name, header, mask=maskM)


def run_stacking_master(fitspath, name,grid_data_file):
    """
    Purpose
    Run function for above function if mask not used
    """
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wavemaster = header['CRVAL1'] + header['CDELT1'] * np.arange(header['NAXIS1'])
    Master_Stacking(fitspath, dataset, wavemaster, grid_data_file, image2DM, name, header, mask=None)
    