
###THIS BINS DATA FOR ALL BINNING METHODS AND IS USED IN GENERAL FUNCTION###


#Currently binning with bin=0.25
#Numpy binary file or FITS binary table: read up on both and attempt both
#Change name of pdf, name of output file, bin size

###Create r 23 grid to have nxm dimensions and update each of the corresponding values, like a table
### m will be 2 and



import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.table import vstack 

#fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
#fitspath='/astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/' 
#pdf_pages = PdfPages(fitspath+'Apr22_R23_O32_bin025_scatter_and_hexbin_MasterGrid.pdf') #open pdf document
#pdf_pages = PdfPages('astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/R23_O32_bin025_scatter_and_hexbin.pdf'
#marker0 = ['b','g','r','m']

#Creates a fits table 
'''for ii in range(1,5):
    file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
    data  = Table(fits.getdata(file1))
    if ii == 1:
        data0 = data
    else:
        data0 = vstack([data0, data])

print 'data0 : ', len(data0)


#For Loop that imports data 

O2 = data0['OII_FLUX_MOD']
O3 = 1.33*data0['OIIIR_FLUX_MOD']
Hb = data0['HB_FLUX_MOD']
R23 = (O2+O3)/Hb
O32 = O3/O2
SNR2 = data0['OII_SNR']
SNR3 = data0['OIIIR_SNR']
SNRH = data0['HB_SNR']

#SNR code: This rules out major outliers by only using specified data
det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                (O2 > 0) & (O3 > 0) & (Hb>0))[0]
'''


def single_grid_O32(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin):
    
    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for O32 

    sort0 = np.argsort(O32)
    y_sort0 = O32[sort0]

    
    #200 galaxies per bin
    n_bins = np.int(len(O32)/galinbin)
    print(n_bins)
    #print(len(y_sort0))

    #Initializing Arrays for outfile later
    N_arr0 = np.zeros((1, n_bins))
    T_arr  = np.zeros((1, n_bins), dtype=object)
    

    #Bin starts and stops initializing
    bin_start = np.zeros(n_bins)
    bin_end   = np.zeros(n_bins)

    #Sets the bins 
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*galinbin]
        bin_end[ii]   = y_sort0[(ii+1)*galinbin-1]
        if ii == max(n_bins):  bin_end[ii] = max(y_sort0)
        #print bin_start[ii] , bin_end[ii]

    #Organizes data into bins
    for oo in range(n_bins):
        idx_arr = np.where((O32>= bin_start[oo]) & (O32<= bin_end[oo]))[0]
        N_arr0[0,oo] += len(idx_arr) 
        T_arr[0,oo]   = idx_arr
    #print'N_arr0:', N_arr0
    #print 'T_arr:', T_arr

    #Plotting
    fig, ax = plt.subplots()
    x = R23
    y = O32
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x1 = x[finite0]
    y1 = y[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    hlines = np.log10(bin_start)
    hlines_end = np.log10(bin_end)
    #print 'start:', hlines
    #print 'end:', hlines_end
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for pp in range(len(bin_start)): plt.axhline(y =hlines[pp] , linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end)): plt.axhline(y =hlines_end[ll], linewidth= 0.3, color= 'g')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    R23_grid = [np.min(R23)]
    
    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid , O32_grid=bin_start, N_arr0=N_arr0)

    fig.clear()

def single_grid_R23(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, galinbin, adaptive='False'):   
    #O2_det3, O3_det3, Hb_det3,

    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    sort0 = np.argsort(R23)
    R23_sort0 = R23[sort0]

    
    #200 galaxies per bin
    #n_bins_old = np.int(len(R23)/galinbin)
    n_bins = len(galinbin)
    print(n_bins)

    #Initializing Arrays for outfile later
    #O32_grid    = np.zeros((n_bins,1))
    R23_grid  = np.zeros((n_bins,1))
    N_arr0 = np.zeros((n_bins,1))
    T_arr  = np.zeros((n_bins,1), dtype = object)
    

    #Bin starts and stops initializing
    bin_start = np.zeros(n_bins, dtype =np.int)
    bin_end   = np.zeros(n_bins, dtype =np.int)

    for ii in range(n_bins):
        #print(2*n_bins)
        if ii == 0:
            bin_start[ii] = 0
        else:
            bin_start[ii]= bin_end[ii-1]+1
        if ii == n_bins-1:
            bin_end[ii] = len(R23_sort0)-1
        else:
            print(galinbin[ii]+bin_start[ii])
            bin_end[ii] = galinbin[ii]+bin_start[ii]-1   
        print('Bin Start:', bin_start[ii] , 'Bin end:', bin_end[ii])

        idx_arr = np.where((R23>= R23_sort0[bin_start[ii]]) & (R23<= R23_sort0[bin_end[ii]]))[0]
        R23_grid[ii,0] = R23_sort0[bin_start[ii]]
        N_arr0[ii,0] += len(idx_arr)
        T_arr[ii,0]   = idx_arr
        print('Start:', R23_grid[ii])   ###This could be incorrect so I need to check this again###
    '''
    #Sets the bins with R23_sorting
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*galinbin]
        bin_end[ii]   = y_sort0[(ii+1)*galinbin-1]
        if ii == (n_bins-1):  bin_end[ii] = np.max(y_sort0)-1
        print bin_start[ii] , bin_end[ii]'''

    #Set bins with indexing
    

    #Organizes data into bins

        
    #print N_arr0

    #Plotting
    fig, ax = plt.subplots()
    x = R23
    y = O32
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x1 = x[finite0]
    y1 = y[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    hlines = np.log10(bin_start)
    hlines_end = np.log10(bin_end)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    for pp in range(len(bin_start)): plt.axvline(x = hlines[pp], linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end)): plt.axvline(x =hlines_end[ll], linewidth= 0.3, color= 'g')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    O32_grid = [np.min(O32)]
    print('R23_grid:', bin_start)
    print('O32_grid:', O32_grid)
    print('T_arr:', T_arr)
    print('N_arr0:', N_arr0)

    np.savez(outfile, T_arr=T_arr, R23_grid=bin_start, O32_grid=O32_grid, N_arr0=N_arr0)

    fig.clear()

def making_Grid(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, R23_bin, O32_bin):

    fig1, ax1 = plt.subplots() 

    xlim = [0.4,50]
    ylim = [0.1,20]

    
    R23_bin  = R23_bin
    O32_bin  = O32_bin
    R23_grid = np.arange(np.log10(xlim[0]), np.log10(xlim[1])+R23_bin, R23_bin)
    O32_grid = np.arange(np.log10(ylim[0]), np.log10(ylim[1])+R23_bin, O32_bin)
    N_arr0 = np.zeros( (len(R23_grid),len(O32_grid)), dtype=np.float)

    N = len(R23_grid)
    M = len(O32_grid)
    T_arr = np.zeros((N, M), dtype = object)

    #Plotting
    label0 = 'Field = data0[det3], N='+str(len(det3))
    x = np.log10(R23)
    y = np.log10(O32)
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x = x[finite0]
    y = y[finite0]
    scatter = ax1.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1, label=label0)
    x0 = x.tolist()
    y0 = y.tolist()

    for jj in range(len(R23_grid)):
        for kk in range(len(O32_grid)):
            array= np.where((x < R23_grid[jj]+R23_bin) & (x >= R23_grid[jj]) &
                            (y < O32_grid[kk]+O32_bin) & (y >= O32_grid[kk]))[0]
            print(array)
            N_arr0[jj,kk]    += len(array)
            T_arr[jj,kk] = array

    ax1.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax1.set_xlabel(r'log($R_{23}$)')
    ax1.set_ylabel(r'log($O_{32}$)')
    ax1.set_xlim(np.log10(xlim))
    ax1.set_ylim(np.log10(ylim))
    ax1.minorticks_on()
    ax1.legend(loc='upper left', numpoints=3)     
    fig1.set_size_inches(9,9)

    pdf_pages.savefig(fig1)

    fig2 = plt.figure()
    ax2 = plt.gca() # ax = plt.subplots() #plt.gcf()
    cm= plt.cm.get_cmap('Blues')

    #Colorbar and hexbin plotting

    tabmastergrid = Table([x0,y0])
    asc.write(tabmastergrid, fitspath+'testmastergrid.tbl', format='fixed_width_two_line')
    hb0 = ax2.hexbin(x0, y0, gridsize=(len(R23_grid),len(O32_grid)), cmap= cm)
    cbaxes= fig2.add_axes([0.135,0.20, 0.75, 0.025])
    cb= fig2.colorbar(hb0, cax=cbaxes, ticks=[0.,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300], orientation= 'horizontal')
    cb.set_label('density')
    #fig.tight_layout()

    ax2.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1, label=label0)
    ax2.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax2.set_xlabel(r'log($R_{23}$)')
    ax2.set_ylabel(r'log($O_{32}$)')
    ax2.set_xlim(np.log10(xlim))
    ax2.set_ylim(np.log10(ylim))
    ax2.minorticks_on()
    ax2.legend(loc='upper left', numpoints=3)

    fig2.set_size_inches(8,8)
    pdf_pages.savefig(fig2)

    pdf_pages.close()

    fig1.clear()
    fig2.clear()

    #outfile = 'Arrays_R23O32bin025MasterGrid.npz' 
    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid, N_arr0=N_arr0)


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






def n_times_binned(fitspath, pdf_pages, outfile, n_split, individual_ID, R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive=False):
    dataset = 'n_Bins'
    ##R23 and O32 are going to be log values --> check to make sure that this doesn't affect the binning
    #One_dimensional binning for R23 followed by each bin being split in O32 in n_split bins 
    '''increase the number of galaxies as R23 increases'''
    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    #R23 = 10**(logR23)
    #O32 = 10**(logO32)
    sortR23 = np.argsort(R23)
    R23_sort0 = R23[sortR23]


    n_bins = len(galinbin)
    
    print(n_bins, 'n_split',n_split)
    n_bins_range = np.arange(0,n_split*n_bins,1)


    #Initializing Arrays for Grid stacking
    Number_inbin = np.zeros((n_bins,n_split))   #Used to be N_arr0
    locator  = np.zeros((n_bins,n_split), dtype = object) #Used to be T_arr
    O32_minimum    = np.zeros((n_bins,n_split))    #Used to be O32_grid
    R23_minimum    = np.zeros((n_bins,n_split))      #Used to be R23_grid
    O32_median = np.zeros((n_bins,n_split))
    R23_median = np.zeros((n_bins,n_split))
    O32_max    = np.zeros((n_bins,n_split))
    R23_max    = np.zeros((n_bins,n_split))
    xBar  = np.zeros(n_split*n_bins)
    yBar  = np.zeros(n_split*n_bins)
    area  = np.zeros(n_split*n_bins)
    Bin_number = np.zeros(len(data3), dtype = int) #Used to be N_bin

    R23_minall = np.zeros(len(R23))
    O32_minall = np.zeros(len(R23))
    R23_avgall = np.zeros(len(R23))
    O32_avgall = np.zeros(len(R23))
    R23_medall = np.zeros(len(R23))
    O32_medall = np.zeros(len(R23))
    R23_maxall = np.zeros(len(R23))
    O32_maxall = np.zeros(len(R23))

    

    #Bin starts and stops initializing
    bin_start = np.zeros(n_bins, dtype=np.int)
    bin_end   = np.zeros(n_bins, dtype=np.int)

    print("bin start:", bin_start)
    print("bin end:" , bin_end)
    
    print("galinbin:", galinbin)
    
    for ii in range(n_bins):
        print('ii: ', ii) 
        if ii == 0:
            bin_start[ii] = 0
        else:
            bin_start[ii]= bin_end[ii-1]+1
        if ii == n_bins-1:
            bin_end[ii] = len(R23_sort0)-1
        else:
            #print galinbin[ii]+bin_start_1[ii]
            bin_end[ii] = galinbin[ii]+bin_start[ii]-1   
        print('Bin Start:', bin_start[ii] , 'Bin end:', bin_end[ii])


        
        R23_idx = np.where((R23>= R23_sort0[bin_start[ii]]) & (R23<= R23_sort0[bin_end[ii]]))[0]
        # There is probably a more streamline way to do this but let's start with version 1 first then update later
        # @{

        R23_inbins = R23[R23_idx]
        O32_inbins = O32[R23_idx]          #O32 relative to the R23_index/ O32 sorted into the R23 bins
        O32_index = np.argsort(O32_inbins) #Sort the O32 in their bins so that we have their locations
        sortO32 = O32_inbins[O32_index]    #Take the O32 in each bin and organizes them based on the index in the previous line
        
        
        
        # }

        
        # The following lines could go into a defintion called "evenSplit( ... )" therefore
        # allowing you to create a different splitting method definition called "optimalSplit( ... )"
        # or something to that effect if so desired.
        # @{
        n_subbins = np.int(np.floor(float(len(sortO32))/n_split))
        subbin_arr = np.ones(n_split, dtype = int)*n_subbins #We have to take the length of the range of n_split because n_split is of type int
        n_remainder = len(sortO32) - n_subbins*n_split
        backwardsIdx = -1
        forwardsIdx = 0
        
        # Sorting all the remainders into the correct bins: Use backwardsIdx or forwardsIdx to start adding spectra
        #from the front or the back of the list of bins 
        for rr in range(n_remainder):
            subbin_arr[backwardsIdx] = n_subbins+1
            backwardsIdx -= 1
        subbin_arr[backwardsIdx] = n_subbins
        # }
        
        startIdx = 0
        endIdx = subbin_arr[0]
        for jj in range(n_split):
            # Let's grab all O32 values
            O32_values_perbin = sortO32[startIdx:endIdx]
            print(O32_values_perbin)
            O32_inbins_idx = O32_index[startIdx:endIdx]

            #This index gives the positions of all the R23 values relative 
            #to the O32 values over the entire R23 bin
            N_bin_idx = R23_idx[O32_inbins_idx]        

            #Gives the bin number for each spectra
            Bin_number[N_bin_idx] = (ii*n_split)+jj  

            # Now let's start storting our data into variables to use later

            #First two map minimum R23 and O32 measure to all individual spectra in bin 
            #Second two values give the lowest O32 and R23 value set for each bin
            R23_minall[N_bin_idx]= R23_sort0[bin_start[ii]]
            O32_minall[N_bin_idx]= O32_values_perbin[0]
            R23_minimum[ii,jj]   = R23_sort0[bin_start[ii]]
            O32_minimum[ii,jj]   = O32_values_perbin[0]

            #First two map median R23 and O32 measure to all individual spectra in bin 
            #Second two give the median R23 and O32 value for each bin
            R23_medall[N_bin_idx]= np.median(R23[N_bin_idx])
            O32_medall[N_bin_idx]= np.median(O32_values_perbin)      
            R23_median[ii,jj]    = np.median(R23[N_bin_idx])
            O32_median[ii,jj]    = np.median(O32_values_perbin)     

            #First two map maximum R23 and O32 measure to all individual spectra in bin 
            #Second two give the maximum R23 and O32 value for each bin
            R23_maxall[N_bin_idx]= np.max(R23[N_bin_idx])
            O32_maxall[N_bin_idx]= np.max(O32[N_bin_idx])
            R23_max[ii,jj]       = np.max(R23[N_bin_idx])
            O32_max[ii,jj]       = np.max(O32[N_bin_idx])
            
            
            #Gives the number of galaxies in each bin
            Number_inbin[ii,jj]  += len(O32_values_perbin)                      

            #Gives the index (location numbers) for the spectra in each bin
            #and is used later loop over and get the galaxies 
            locator[ii,jj]   = N_bin_idx                                    

            #Maps average R23 measure to all individual spectra in bin 
            R23_avgall[N_bin_idx]= np.average(R23[N_bin_idx])
            #Gives the average R23 value for each bin
            xBar[(ii*n_split)+jj]= np.average(R23[N_bin_idx])
            
              
            #Maps average O32 measure to all individual spectra in bin
            O32_avgall[N_bin_idx]= np.average(O32_values_perbin)
            #Gives the average O32 value for each bin
            yBar[(ii*n_split)+jj]= np.average(O32_values_perbin)

            #Gives the number of galaxies in each bin
            area[(ii*n_split)+jj]= len(O32_values_perbin)                 

            # Now we can shift our window over for the next split like the counting index the Stacking code
            startIdx = endIdx
            endIdx = startIdx+subbin_arr[jj]+1


    O32_lowlimit = O32_minimum.reshape(n_split*n_bins) #Used to be O32_grids
    R23_lowlimit = R23_minimum.reshape(n_split*n_bins) #Used to be R23_grids
    O32_medians = O32_median.reshape(n_split*n_bins)
    R23_medians = R23_median.reshape(n_split*n_bins)
    R23_maxval = R23_max.reshape(n_split*n_bins)
    O32_maxval = O32_max.reshape(n_split*n_bins)
    
        
    #Plotting
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    vlines = np.log10(R23_lowlimit)
    hlines = np.log10(O32_lowlimit)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    #for pp in range(n_split*len(bin_start)): plt.axvline(x = vlines[pp], linewidth= 0.3, color= 'k')

    #print "R23_grids:", np.log10(R23_lowlimit)
    #print "O32_grids:", np.log10(O32_lowlimit)

    for jj in range(len(O32_lowlimit)):
        xmin = vlines[jj]
        if jj <= (len(O32_lowlimit)-n_split-1): xmax = vlines[jj+n_split]
        else: xmax = np.log10(max(R23))
        #print "jj, xmin, xmax, hlines[jj], vlines[jj]", jj, xmin, xmax, hlines[jj], vlines[jj]
        plt.axvline(x = vlines[jj], linewidth= 0.3, color= 'k')
        
        x_value = [xmin,xmax]
        y_value = [hlines[jj], hlines[jj]]
        y_average = [yBar[jj],yBar[jj]]
        plt.plot(x_value,y_value, linewidth= 0.3, color= 'b')
        plt.plot(x_value, y_average, linewidth= 0.3, color= 'g')



        #plt.xlim(-0.3, 1)
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()


    np.savez(outfile, locator=locator, R23_minimum=R23_minimum, O32_minimum=O32_minimum, Number_inbin=Number_inbin)

    n1 = ('bin_ID' , 'logR23_min', 'logO32_min', 'logR23_avg','logO32_avg',
          'logR23_median','logO32_median','logR23_max','logO32_max', 'N_stack')
    tab1 = Table([n_bins_range, R23_lowlimit, O32_lowlimit,xBar, yBar,
                  R23_medians, O32_medians, R23_maxval, O32_maxval, area], names = n1)
    asc.write(tab1, fitspath+'/bin_info.tbl', format='fixed_width_two_line')   #used to be called +dataset+'_binning_averages.tbl
    
    fig.clear()

    n2=('logR23', 'logO32', 'OIII_5007_S/N', 'bin_ID','ID',
        'logR23_min', 'logO32_min', 'logR23_avg','logO32_avg',
        'logR23_median','logO32_median','logR23_max','logO32_max')
    tab2= Table([R23, O32, SNR3, Bin_number, individual_ID,
                 R23_minall,O32_minall,R23_avgall,O32_avgall,
                 R23_medall,O32_medall,R23_maxall, O32_maxall], names=n2)
    asc.write(tab2, fitspath+'/individual_bin_info.tbl', format='fixed_width_two_line') #used to be + dataset+'_2d_binning_datadet3.tbl

    '''n3 = ('ID' , 'R23_grid', 'O32_grid')
    tab1 = Table([n_bins_range, R23_grid, O32_grid], names = n3)
    asc.write(tab1, fitspath+'/Double_Bin_grid_values.tbl', format='fixed_width_two_line')'''
###Create another ascii table with the R23_grid and O32_grid values for plots 
