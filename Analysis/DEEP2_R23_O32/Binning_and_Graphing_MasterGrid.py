#Currently binning with bin=0.25
#Numpy binary file or FITS binary table: read up on both and attempt both
#Change name of pdf, name of output file, bin size

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


def single_grid_O32(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3,galinbin):
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

def single_grid_R23(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3,galinbin):


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
        if ii == max(n_bins):  bin_end[ii] = max(y_sort0)
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

def making_Grid(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3, R23_bin, O32_bin):

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


def two_times_binned(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3,galinbin):

    #One_dimensional binning for R23 followed by each bin being split in half one with high O32 and one with low O32

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
    bin_start_1 = np.zeros(n_bins)
    bin_end_1   = np.zeros(n_bins)
    bin_start_2 = np.zeros(n_bins)
    bin_end_2   = np.zeros(n_bins)

    #Sets the bins 
    for ii in range(n_bins):
        bin_start_1[ii] = y_sort0[ii*galinbin]
        bin_end_1[ii]   = y_sort0[(ii+1)*galinbin-1]
        print bin_start_1[ii] , bin_end_1[ii]


    #Splitting bins in half by O32 value
    sort_O32 = np.argsort(O32)
    x_sort = O32[sort_O32]
    for aa in range(n_bins):
        bin_start_2[ii] = x_sort[ii*galinbin]
        bin_end_2[ii]   = x_sort[(ii+1)*galinbin-1]
        print bin_start_2[ii] , bin_end_2[ii]
        
    #Organizes data into bins
    for oo in range(n_bins):
        for kk in range(n_bins):
            idx_arr = np.where((R23>= bin_start_1[oo]) & (R23<= bin_end_1[oo])
                           & (O32>= bin_start_2[kk]) & (O32<= bin_end_2[kk]))[0]
            print idx_arr
            N_arr0[oo,kk] += len(idx_arr)
            T_arr[oo,kk]   = idx_arr
            print 'N_arr0:', N_arr0
            print 'T_arr:', T_arr

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
    vlines = np.log10(bin_start_2)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    for pp in range(len(bin_start_1)): plt.axvline(x = hlines[pp], linewidth= 0.3, color= 'k')
    for ll in range(len(bin_end)): plt.axvline(x =hlines_end[ll], linewidth= 0.3, color= 'g')
    for jj in range(len(bin_start_2)): plt.axhline(y =vlines[jj], linewidth= 0.3, color= 'b')
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    

    np.savez(outfile, T_arr=T_arr, R23_grid=bin_start_1, O32_grid=bin_start_2, N_arr0=N_arr0)

    fig.clear()
