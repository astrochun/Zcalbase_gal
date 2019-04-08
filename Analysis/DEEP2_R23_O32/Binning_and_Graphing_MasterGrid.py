
###THIS BINS DATA FOR ALL BINNING METHODS AND IS USED IN GENERAL FUNCTION###


#Currently binning with bin=0.25
#Numpy binary file or FITS binary table: read up on both and attempt both
#Change name of pdf, name of output file, bin size

###Create r 23 grid to have nxm dimensions and update each of the corresponding values, like a table
### m will be 2 and
'''R23 grid np.zeros ((5,2)) for rr in range(R23_grid.shape[0])) print rr
for oo in range ((R23_grid.shape[1])) print oo
define n_N = R23_grid.shape[0] or [1]

R23_grid[0][:] = 3.0'''


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
    #O2_det3, O3_det3, Hb_det3
    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for O32 

    sort0 = np.argsort(O32)
    y_sort0 = O32[sort0]

    
    #200 galaxies per bin
    n_bins = np.int(len(O32)/galinbin)
    print n_bins
    #print len(y_sort0)

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
        print bin_start[ii] , bin_end[ii]

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
    for pp in range(len(bin_start)): plt.axhline(y =hlines[pp] , linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end)): plt.axhline(y =hlines_end[ll], linewidth= 0.3, color= 'g')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    R23_grid = [np.min(R23)]
    
    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid , O32_grid=bin_start, N_arr0=N_arr0)

    fig.clear()

def single_grid_R23(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, galinbin):   
    #O2_det3, O3_det3, Hb_det3,

    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    sort0 = np.argsort(R23)
    y_sort0 = R23[sort0]

    
    #200 galaxies per bin
    n_bins = np.int(len(R23)/galinbin)
    print n_bins

    #Initializing Arrays for outfile later
    N_arr0 = np.zeros((n_bins,1))
    T_arr  = np.zeros((n_bins,1), dtype = object)
    

    #Bin starts and stops initializing
    bin_start = np.zeros(n_bins)
    bin_end   = np.zeros(n_bins)

    #Sets the bins 
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*galinbin]
        bin_end[ii]   = y_sort0[(ii+1)*galinbin-1]
        if ii == (n_bins-1):  bin_end[ii] = np.max(y_sort0)-1
        print bin_start[ii] , bin_end[ii]

    #Organizes data into bins
    for oo in range(n_bins):
        idx_arr = np.where((R23>= bin_start[oo]) & (R23<= bin_end[oo]))[0]
        N_arr0[oo,0] += len(idx_arr)
        T_arr[oo,0]   = idx_arr
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
    for pp in range(len(bin_start)): plt.axvline(x = hlines[pp], linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end)): plt.axvline(x =hlines_end[ll], linewidth= 0.3, color= 'g')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    O32_grid = [np.min(O32)]

    np.savez(outfile, T_arr=T_arr, R23_grid=bin_start, O32_grid=O32_grid, N_arr0=N_arr0)

    fig.clear()

def making_Grid(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, R23_bin, O32_bin):
    #O2_det3, O3_det3, Hb_det3,
    fig1, ax1 = plt.subplots() #plt.gcf()

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
            print array
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


def two_times_binned(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin):
    #O2_det3, O3_det3, Hb_det3
    #One_dimensional binning for R23 followed by each bin being split in half one with high O32 and one with low O32

    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    sort0 = np.argsort(R23)
    y_sort0 = R23[sort0]

    
    #200 galaxies per bin
    #
    #for oo in range(len(galinbin)):
    n_bins = np.int(len(R23)/galinbin)
    print n_bins
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
    bin_start_1 = np.zeros(n_bins)
    bin_end_1   = np.zeros(n_bins)
    

    #Sets the bins
    #galinbin = [100,150, 250, 300, 300, 200, 205, 102, 102]
    '''increase the number of galaxies as R23 increases'''
    #if len(galinbin) == n_bins:
        #print 'Length of galinbin list equal to number of calculated bins'

    for ii in range(n_bins):
        print galinbin

        '''if ii == 0: bin_start_1 = 0
        else: bin_start_1[ii]= y_sort0[ii*galinbin[ii-1]]
        bin_end_1[ii]= y_sort0[galinbin[ii]+bin_start_1]'''

        '''bin_start_1[ii] = y_sort0[ii*galinbin]  #[ii]
        bin_end_1[ii]   = y_sort0[(ii+1)*galinbin-1]  #[ii]'''
        if ii == n_bins-1:  bin_end_1[ii] = np.max(y_sort0)
        #print 'Bin Start:', bin_start_1[ii] , 'Bin end:', bin_end_1[ii]
        #print 'Bin Start:', bin_start_1 , 'Bin end:', bin_end_1
        idx1 = np.where((R23>= bin_start_1[ii]) & (R23<= bin_end_1[ii]))[0]
        #idx1 = np.where((R23>= bin_start_1) & (R23<= bin_end_1))[0]
        med_val = np.median(O32[idx1])

        idx2 = np.where((R23>= bin_start_1[ii]) & (R23<= bin_end_1[ii]) & (O32 <= med_val))[0]
        idx3 = np.where((R23>= bin_start_1[ii]) & (R23<= bin_end_1[ii]) & (O32 > med_val))[0]

        '''idx2 = np.where((R23>= bin_start_1) & (R23<= bin_end_1) & (O32 <= med_val))[0]
        idx3 = np.where((R23>= bin_start_1) & (R23<= bin_end_1) & (O32 > med_val))[0]'''

            
        '''O32_grid[ii*2]   = np.median(O32[idx2])    
        O32_grid[ii*2+1] = np.median(O32[idx3])
        R23_grid[ii*2] = bin_start_1[ii]
        R23_grid[ii*2+1] = bin_start_1[ii]'''
        
        O32_grid[ii,0]   = np.median(O32[idx2])    
        O32_grid[ii,1] = np.median(O32[idx3])

        R23_grid[ii,0] = bin_start_1[ii]
        R23_grid[ii,1] = bin_start_1[ii]
        
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

        O32_values[ii*2]   = np.median(O32[idx2])    
        O32_values[ii*2+1] = np.median(O32[idx3])
        R23_values[ii*2] = bin_start_1[ii]
        R23_values[ii*2+1] = bin_start_1[ii]
                
        '''if ii== n_bins-2 or ii ==n_bins-1: 
            O32_grid[ii,0]=(np.median(O32[idx2]))/2
        O32_grid[ii,1]=(np.median(O32[idx2]))/2
        O32_grid[ii,2]=(np.median(O32[idx2]))/2      but obviously this we would need to redefine our initialization of areas'''
                

        
    #Plotting
    fig, ax = plt.subplots()
    x = R23
    y = O32
    finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    x1 = x[finite0]
    y1 = y[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    hlines = np.log10(bin_start_1)
    hlines_end = np.log10(bin_end_1)
    #vlines = np.log10(bin_start_2)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    for pp in range(len(bin_start_1)): plt.axvline(x = hlines[pp], linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end_1)): plt.axvline(x =hlines_end[ll], linewidth= 0.3, color= 'g')
    #for jj in range(len(bin_start_2)): plt.axhline(y =vlines[jj], linewidth= 0.3, color= 'b')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    
    '''print 'bin:', len(n_bins_range)
    print 'R23:', len(R23_grid), R23_grid
    print 'bin_start_1', len(bin_start_1), bin_start_1
    print 'O32:', len(O32_grid)'''

    np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid, N_arr0=N_arr0)

    n1 = ('ID' , 'R23_value', 'O32_value', 'xBar','yBar', 'area')
    tab1 = Table([n_bins_range, R23_values, O32_values,xBar, yBar,area], names = n1)
    asc.write(tab1, fitspath+'/Double_Bin_binning_averages.tbl', format='fixed_width_two_line')
    
    fig.clear()

    n2=('R23', 'O32', 'N_bin')
    tab2= Table([R23, O32, N_bin], names=n2)
    asc.write(tab2, fitspath+'/Double_Bin_2d_binning_datadet3.tbl', format='fixed_width_two_line')

###Create another ascii table with the R23_grid and O32_grid values for plots 
