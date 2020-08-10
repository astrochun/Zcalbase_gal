
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table

def two_times_binned(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive=False):
    print(fitspath)
    #O2_det3, O3_det3, Hb_det3
    #One_dimensional binning for R23 followed by each bin being split in half one with high O32 and one with low O32

    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    sort0 = np.argsort(R23)
    R23_sort0 = R23[sort0]

    print("Running :)")
    #200 galaxies per bin
    #
    #for oo in range(len(galinbin)):
    n_bins = len(galinbin)
    '''if adaptive == True: n_bins = len(galinbin)
    if adaptive == False: n_bins = np.int(len(R23)/galinbin)'''
    
    print(n_bins)  
    n_bins_range = np.arange(0,2*n_bins,1)

    #Initializing Arrays for outfile later if using Voronoi Stacking
    '''N_arr0 = np.zeros(2*n_bins)
    T_arr  = np.zeros((n_bins*2), dtype = object)
    O32_grid    = np.zeros(2*n_bins)
    R23_grid  = np.zeros(2*n_bins)
    xBar  = np.zeros(2*n_bins)
    yBar  = np.zeros(2*n_bins)
    area  = np.zeros(2*n_bins)
    N_bin = np.zeros(len(data3), dtype = int)'''

    #Initializing Arrays for Grid stacking
    N_arr0 = np.zeros((n_bins,2))
    T_arr  = np.zeros((n_bins,2), dtype = object)
    O32_grid    = np.zeros((n_bins,2))
    R23_grid  = np.zeros((n_bins,2))
    xBar  = np.zeros(2*n_bins)
    yBar  = np.zeros(2*n_bins)
    area  = np.zeros(2*n_bins)
    N_bin = np.zeros(len(data3), dtype = int)

    O32_values    = np.zeros(2*n_bins)
    R23_values  = np.zeros(2*n_bins)
    

    #Bin starts and stops initializing
    bin_start_1 = np.zeros(n_bins, dtype=np.int)
    bin_end_1   = np.zeros(n_bins, dtype=np.int)

    #Sets the bins
    #galinbin = [100,150, 250, 300, 300, 200, 205, 102, 102]
    '''increase the number of galaxies as R23 increases'''
    #if len(galinbin) == n_bins:
        #print('Length of galinbin list equal to number of calculated bins')

    print(galinbin)
    for ii in range(n_bins):
        #print(2*n_bins)
        if ii == 0:
            bin_start_1[ii] = 0
        else:
            bin_start_1[ii]= bin_end_1[ii-1]+1
        if ii == n_bins-1:
            bin_end_1[ii] = len(R23_sort0)-1
        else:
            print(galinbin[ii]+bin_start_1[ii])
            bin_end_1[ii] = galinbin[ii]+bin_start_1[ii]-1   
        print('Bin Start:', bin_start_1[ii] , 'Bin end:', bin_end_1[ii])


        '''bin_start_1[ii] = y_sort0[ii*galinbin]  #[ii]
        bin_end_1[ii]   = y_sort0[(ii+1)*galinbin-1]  #[ii]'''


        #if ii == n_bins-1:  bin_end_1[ii] = np.max(R23_sort0)
        #print('Bin Start:', bin_start_1[ii] , 'Bin end:', bin_end_1[ii])

        idx1 = np.where((R23>= R23_sort0[bin_start_1[ii]]) & (R23<= R23_sort0[bin_end_1[ii]]))[0]
        #idx1 = np.where((R23>= bin_start_1) & (R23<= bin_end_1))[0]
        med_val = np.median(O32[idx1])

        idx2 = np.where((R23>= R23_sort0[bin_start_1[ii]]) & (R23<= R23_sort0[bin_end_1[ii]]) & (O32 <= med_val))[0]
        idx3 = np.where((R23>= R23_sort0[bin_start_1[ii]]) & (R23<= R23_sort0[bin_end_1[ii]]) & (O32 > med_val))[0]

        '''idx2 = np.where((R23>= bin_start_1) & (R23<= bin_end_1) & (O32 <= med_val))[0]
        idx3 = np.where((R23>= bin_start_1) & (R23<= bin_end_1) & (O32 > med_val))[0]'''

            
        O32_grid[ii,0]   = np.median(O32[idx2])    
        O32_grid[ii,1] = np.median(O32[idx3])
        R23_grid[ii,0] = R23_sort0[bin_start_1[ii]]
        R23_grid[ii,1] = R23_sort0[bin_start_1[ii]]
        
        '''O32_grid[ii,0]   = np.median(O32[idx2])    
        O32_grid[ii,1] = np.median(O32[idx3])

        R23_grid[ii,0] = bin_start_1[ii]
        R23_grid[ii,1] = bin_start_1[ii]'''
        
        '''N_arr0[ii*2]  += len(idx2)
        N_arr0[ii*2+1]+= len(idx3)
        
        T_arr[ii*2]   = idx2
        T_arr[ii*2+1] = idx3'''

        N_arr0[ii,0]  += len(idx2)
        N_arr0[ii,1]+= len(idx3)
        
        T_arr[ii,0]   = idx2
        T_arr[ii,1] = idx3
        
        N_bin[idx2]= ii*2
        N_bin[idx3]= ii*2+1

        

        xBar[ii*2]= np.log10(np.average(R23[idx2]))
        xBar[ii*2+1]= np.log10(np.average(R23[idx3]))

        yBar[ii*2]= np.log10(np.average(O32[idx2]))
        yBar[ii*2+1]= np.log10(np.average(O32[idx3]))
        
        area[ii*2]= len(idx2)
        area[ii*2+1]= len(idx3)

        
        
                
        '''if ii== n_bins-2 or ii ==n_bins-1: 
            O32_grid[ii,0]=(np.median(O32[idx2]))/2
        O32_grid[ii,1]=(np.median(O32[idx2]))/2
        O32_grid[ii,2]=(np.median(O32[idx2]))/2      but obviously this we would need to redefine our initialization of areas'''
                

    O32_values   = O32_grid.reshape(2*n_bins) 
    R23_values = R23_grid.reshape(2*n_bins)
    print('O32_values', O32_values)
    print('R23_values', R23_values)
        
    #Plotting
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    hlines = np.log10(R23_values)
    #hlines_end = np.log10()
    #vlines = np.log10(bin_start_2)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for pp in range(2*len(bin_start_1)): plt.axvline(x = hlines[pp], linewidth= 0.3, color= 'k')
    #for ll in range(len(bin_end_1)): plt.axvline(x =hlines_end[ll], linewidth= 0.3, color= 'g')
    #for jj in range(len(bin_start_2)): plt.axhline(y =vlines[jj], linewidth= 0.3, color= 'b')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

   
    print("Please work :)")
    
    '''print 'bin:', len(n_bins_range)
    print 'R23:', len(R23_grid), R23_grid
    print 'bin_start_1', len(bin_start_1), bin_start_1
    print 'O32:', len(O32_grid)'''

    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid, N_arr0=N_arr0)

    n1 = ('ID' , 'R23_value', 'O32_median', 'xBar','yBar', 'area')
    tab1 = Table([n_bins_range, R23_values, O32_values,xBar, yBar,area], names = n1)
    asc.write(tab1, fitspath+'/Double_Bin_binning_averages.tbl', format='fixed_width_two_line')
    
    fig.clear()

    n2=('R23', 'O32', 'SN_5007', 'N_bin')
    tab2= Table([R23, O32, SNR3, N_bin], names=n2)
    asc.write(tab2, fitspath+'/Double_Bin_2d_binning_datadet3.tbl', format='fixed_width_two_line')

    '''n3 = ('ID' , 'R23_grid', 'O32_grid')
    tab1 = Table([n_bins_range, R23_grid, O32_grid], names = n3)
    asc.write(tab1, fitspath+'/Double_Bin_grid_values.tbl', format='fixed_width_two_line')'''
###Create another ascii table with the R23_grid and O32_grid values for plots 

