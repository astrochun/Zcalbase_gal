
import numpy as np
import matplotlib.pyplot as plt
#import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
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
import scipy.integrate as integ

import general

#fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_104/'



#single_vs_averaged = fitspath +'single_vs_averaged_emissions.pdf'
#A1_attempt2 = fitspath +'A1_attempt2.pdf'
#A2_attempt2 = fitspath +'A2_attempt2.pdf'
#comparing_composite_values = fitspath +'comparing_composite_values.pdf'
'''
#Raw R23 and O32 (xnode is R23 and ynode is O32) these are average values calculated by voronoi code 10
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi.tbl'
data0 = asc.read(raw)
R23_raw_voronoi = data0['xBar']
O32_raw_voronoi = data0['yBar']'''
'''
#Raw R23 and O32 (xnode is R23 and ynode is O32) these are average values calculated by voronoi code 14
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_14.tbl'
data0 = asc.read(raw)
R23_raw_voronoi = data0['xBar']
O32_raw_voronoi = data0['yBar']

#Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Voronoi Outputs10: R_23 and O_32 values for all galaxies with a column specifying bins 
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar21_voronoi_2d_binning_output.txt'
voronoi = np.genfromtxt(outfilevoronoi)

#Voronoi Outputs14: R_23 and O_32 values for all galaxies with a column specifying bins 
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output_14.txt'
voronoi = np.genfromtxt(outfilevoronoi)

#Grid Inputs and Outputs
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_asc_table_Average_R23_O32_Values.tbl' 
data0 = asc.read(raw)
R23_raw_voronoi = data0['R_23_Average']
O32_raw_voronoi = data0['O_32_Average']

spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_combined_flux_table.tbl'
data1 = asc.read(spectral)'''

def Plotting_Data1(fitspath, dataset, combine_flux_ascii, asc_table1):
    
    line_plot = fitspath +dataset+'_line_ratio_plots.pdf'
    
    #combine_flux_ascii = fitspath+dataset#+'_combined_flux_table.tbl'
    print "### combine_flux_ascii : "+combine_flux_ascii 
    fitted_data = asc.read(combine_flux_ascii)

    print "### asc_table1 : "+asc_table1
    raw_data = asc.read(asc_table1)
    
    OII = fitted_data['OII_3727_Flux_Observed']
    OIII4959 = fitted_data['OIII_4958_Flux_Observed']
    OIII5007 = fitted_data['OIII_5007_Flux_Observed']
    H_BETA = fitted_data['HBETA_Flux_Observed']
    binnum = fitted_data['N_Galaxies']
    print 'binnum:', binnum, len(binnum)
    pdf_pages = PdfPages(line_plot)
    nrows = 4
    ncols = 4

    R23_composite = np.zeros(binnum.shape[0])#np.log10((OII + (1.33*OIII5007))/H_BETA)
    O32_composite = np.zeros(binnum.shape[0]) #np.log10((1.33*OIII5007)/OII)
    for ii in range(len(binnum)):
        R23_comp = np.log10((OII[ii] + (1.33*OIII5007[ii]))/H_BETA[ii])
        O32_comp = np.log10((1.33*OIII5007[ii])/OII[ii])
        print R23_comp, O32_comp
        R23_composite[ii]= R23_comp
        O32_composite[ii]= O32_comp
    
    R23_raw = raw_data['xBar']
    O32_raw = raw_data['yBar']
    binnum_raw = raw_data['area']
    print 'binnum_raw', binnum_raw, len(binnum_raw) 
    #for aa in range(len(binnum)+1): print R23_raw[aa], O32_raw[aa]#, binnum_raw[aa]

    if dataset != 'Grid':
        for rr in range(len(binnum)):
            if binnum[rr] == binnum_raw[rr]:
                print 'equal',binnum[rr], binnum_raw[rr]
            #print 'R23 composite vs raw' R23_composite[rr], R23_raw[rr]
            #print 'O32 composite vs raw' O32_composite[rr], O32_raw[rr]
            #else: print binnum,  binnum_raw +'not equal...NOOOOO'
        
    #R23= np.log10(R23_composite/R23_raw)
    #O32= np.log10(O32_composite/O32_raw)
    print 'R23_raw: as calculated by the grid or voronoi code', R23_raw, 'O32_raw: as calculated by the grid or voronoi code', O32_raw
    print 'R23_composite: as calculated from observations', R23_composite, 'O32_composite: as calculated from observations', O32_composite 
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(R23_raw,R23_composite, marker= 'o', facecolor= 'none', edgecolor ='b',label= 'R23 Ratio: Vornoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset+' Raw vs. Composite for R23')
    #for rr in range(len(data1)):
        #if binnum[rr] <= 30: ax_arr.annotate(str(binnum[rr]), (R23_raw_voronoi[rr],R23_composite[rr]), xycoords='data')
    ax_arr.set_xlabel(r'Raw log($R_{23}$)')
    ax_arr.set_ylabel(r'Composite log($R_{23}$)')

    ax_arr.plot([0.0,1.3], [0.0,1.3], 'k-')
    #ax_arr.plot([0.7962914706935305],[0.826262793371], 'ro')
    
    #plt.draw()
    fig.savefig(pdf_pages, format='pdf')

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(O32_raw,O32_composite, marker= 'o', facecolor= 'none', edgecolor ='b', label= 'O32 Ratio: Vornoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset+'Raw vs. Composite for O32')
    #for oo in range(len(data1)):
        #if binnum[oo] <= 30: ax_arr.annotate(str(binnum[oo]), (O32_raw_voronoi[oo],O32_composite[oo]), xycoords='data')
        
    ax_arr.set_xlabel(r'Raw log($O_{32}$)')
    ax_arr.set_ylabel(r'Composite log($O_{32}$)')

    ax_arr.plot([-1,1.2], [-1,1.2], 'k-')
    #ax_arr.plot([0.09102410573242697],[-0.0346399931757], 'ro')
    
    #plt.draw()
    fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()



####The rest of this code is not fully developed###
#RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'

def average_vals_bin2(): 
    '''for ii in range(1,5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

        #print 'data0 : ', len(data0)
    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']'''

    O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3 = general.get_det3()
    #O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3= general.get_det3()

    O2_det3=O2[det3]
    O3_det3=O3[det3]
    Hb_det3=Hb[det3]
    
    #Calculating R23 and O32 for the orginal data
    '''R23_original = np.zero(len(O2_det3))
    O32_original = np.zero(len(O2_det3))
    for ii in len(O2_det2): 
        R23 = (O2_det3[ii]+O3_det3[ii])/Hb_det3[ii]
        O32 = O3_det3[ii]/O2_det3[ii]
        R23_original[ii] = R23
        O32_original[ii] = O32'''
       
    #Composite Spectrum
    lR23 = voronoi[:,2]
    lO32 = voronoi[:,3]
    binnum = voronoi[:,4]
    binnum0= list(set(binnum))
    
    R23_original = []  #np.zeros(len(binnum0), dtype=np.float64)
    O32_original = []  #np.zeros(len(binnum0), dtype=np.float64)
    #print type(avg_l_R23), type(avg_l_O32)
    #print avg_l_R23
    for ii in binnum0:
        idx = np.array([xx for xx in range(len(det3)) if binnum[xx] == ii])
            
        O2_original = np.mean(O2_det3[idx])/1e-17
        O3_original = np.mean(O3_det3[idx])/1e-17
        Hb_original = np.mean(Hb_det3[idx])/1e-17
        #print 'Avg_O3', avg_O3
        list_R23 = np.log10((O2_original +O3_original)/Hb_original)
        list_O32 = np.log10(O3_original/O2_original)
        #print list_R23, type(list_R23), list_O32, type(list_O32)

        R23_original.append(list_R23)
        O32_original.append(list_O32)
    Plotting_Data1(R23_original, O32_original)
 






#The rest of the in this file does different parts of the above two functions. The functions above are the products of what is below. 

def A1_attempt_two_line_plotting():
    #When I loop over the for loop, the data from all previous loops also graphs, how do I fit that? 
    com_R23 = data1['R_23_Average']
    com_O32 = data1['O_32_Average']
    
    binNum = voronoi[:,4]
    
    
    binNum0 = list(set(binNum))
    R23_arr = np.zeros(len(binNum0))
    O32_arr = np.zeros(len(binNum0))

    pdf_pages= PdfPages(A1_attempt2)
    for ii in range(len(binNum0)):
        fig, ax_arr = plt.subplots()
        random_bin = binNum0[ii]
        bin_idx = [xx for xx in range(len(voronoi)) if binNum[xx] == random_bin]

        R23_arr= voronoi[:,1][bin_idx]
        O32_arr= voronoi[:,2][bin_idx]

        c_R23 = com_R23[ii]
        c_O32 = com_O32[ii]
     
        ax_arr.scatter(R23_arr, O32_arr)
        ax_arr.plot(c_R23, c_O32, '*',markersize=12)
        ax_arr.set_xlabel('R23')
        ax_arr.set_ylabel('O32')
        fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()
    
def A2_attempt_two_line_plotting():
    binnum = list(range(0,27))
    pdf_pages= PdfPages(A2_attempt2)
    scalefact = 1e-17
#Stuff for every galaxies
    for ii in range(1,5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

        #print 'data0 : ', len(data0)
    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']
     
    SNR2 = data0['OII_SNR']
    SNR3 = data0['OIIIR_SNR']
    SNRH = data0['HB_SNR']
    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                    (O2 > 0) & (O3 > 0) & (Hb>0))[0]


    data3 = data0[det3]
    O2 = data3['OII_FLUX_MOD']/scalefact
    O3 = 1.33*data3['OIIIR_FLUX_MOD']/scalefact
    Hb = data3['HB_FLUX_MOD']/scalefact

   
    com_O2 = data1['OII_3727_Flux_Observed']
    com_O3 = data1['OIII_4958_Flux_Observed']
    com_Hb = data1['HBETA_Flux_Observed']

    

    '''for ii in binnum:
        fig, ax_arr = plt.subplots()
        ax_arr.scatter(O2,O3) #These numbers are not binned and so I need to figure out how to call the binning 
        ax_arr.plot(com_O2[ii],com_O3[ii], '*', markersize=12)
        fig.savefig(pdf_pages, format='pdf')

    for ii in range(binnum):
        

    pdf_pages.close()'''




def deal_with_later():
#compare the averages for the individual lines that were calculated on Friday with the large table from zm_general
    O2_observed = np.log10(data1['OII_3727_Flux_Observed'])
    O3_observed = np.log10(data1['OIII_5007_Flux_Observed'])
    Hb_observed = np.log10(data1['HBETA_Flux_Observed'])
    #print O2_observed
    binnum = data1['N_Galaxies']
    pdf_pages= PdfPages(single_vs_averaged)
    print avg_O2, avg_O2.size, avg_O3, avg_O3.size

    #O2
    '''row = rr / nrows % ncols
        col = rr % ncols
        #print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]'''
    nrows = 4
    ncols = 4 
    
    fig, ax_arr = plt.subplots()
    
    for ii in range(O2_observed.shape[0]): 
   
      ax_arr.scatter(O2_observed, O3_observed, label= 'O2_observed')
      ax_arr.plot(avg_O2, avg_O3)
      #ax_arr.plot(aver_O2, aver_O3)
      #ax_arr.legend(loc=0)
      #for bb in range(len(data1)):
      #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
      #ax_arr.set_xlabel()
      ax_arr.set_ylabel('O2_Flux')

      #ax_arr.plot([0.0,1.3], [0.0,1.3], 'k-')

      fig.savefig(pdf_pages, format='pdf')

    '''

    #O3
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(binnum,O3_observed, label= 'O3')
    ax_arr.scatter(binnum,avg_O3)
    ax_arr.legend(loc=0)
    #for bb in range(len(data1)):
        #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
    ax_arr.set_xlabel('N_Galaxies')
    ax_arr.set_ylabel('O3_Flux')
    fig.savefig(pdf_pages, format='pdf')

    #Hb
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(binnum,Hb_observed, label= 'Hb')
    ax_arr.scatter(binnum,avg_Hb)
    ax_arr.legend(loc=0)
    #for bb in range(len(data1)):
        #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
    ax_arr.set_xlabel('N_Galaxies')
    ax_arr.set_ylabel('Hb_Flux')
    fig.savefig(pdf_pages, format='pdf')
'''
    pdf_pages.close()







#R23 = ([OII] + [OIII]4959 +[OIII]5007)/H_BETAemission
#O23 = ([OIII]4959 +[OIII]5007)/[OII]
#First you need to find the R23 and O32 values for every line for all emission plots (1-27) using the ascii table. Use the spectral fluxs you calculated
#Then calculate the R23 and O32 values from the raw data(x_node and y_node values yes from Stacking_Voronio) 
#Finally generate plots with the raw R23 and the spectral R23 (same for O32) == you should find a one to one ratio: plot the log vs. log 



#Random Bin xBar and yBar Calculations
#Calculate the average R23 and O32 values for a random bin == import voronoi_2d_binning_output.txt file and for a random bin fidn the number of galaxies and sum??

def checking_R23_O32():
    outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
    voronoi = asc.read(outfilevoronoi)
    binNum = voronoi['col3']

    binNum0 = list(set(binNum))
    R23_arr = np.zeros(len(binNum0))
    O32_arr = np.zeros(len(binNum0))
    
    for ii in range(len(binNum0)):
        random_bin = binNum0[ii]
        bin_idx = [xx for xx in range(len(voronoi)) if binNum[xx] == random_bin]


        R23_arr[ii] = np.average(voronoi['col1'][bin_idx])
        O32_arr[ii] = np.average(voronoi['col2'][bin_idx])
    n= ('R_23_Average', 'O_32_Average', 'bin_number')
    tab1 = Table([R23_arr, O32_arr,binNum0], names=n)
    asc.write(tab1, fitspath+'avg_calcul_values.tbl', format='fixed_width_two_line')


#initialize arrays and populate 

def comparing_tables():
    outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
    voronoi = asc.read(outfilevoronoi)
    avg_calcul_tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/avg_calcul_values.tbl'
    calcul_tab = asc.read(avg_calcul_tab)
    #data1 is the other table
    pdf_pages = PdfPages('/Users/reagenleimbach/Desktop/Zcalbase_gal/comparing_tables.pdf')
    
    binNum = voronoi['col3']
    print len(binNum)
    
    #for ii in range(len(binNum)):
    #computer calculated - checking_R23_O32 calucation 
    diff_R23 = data1['R_23_Average']- calcul_tab['R_23_Average']
    diff_O32 = data1['O_32_Average']- calcul_tab['O_32_Average']
    #if diff_R23 != 0: print diff_R23
    #if diff_O32 != 0: print diff_O32

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(diff_R23,diff_O32)
    #ax_arr.scatter(binNum[ii],diff_O32)
    #ax_arr.set_xlabel(r'Raw log($R_{23}$)')
    #ax_arr.set_ylabel(r'Observed log($R_{23}$)')
    
    '''
    fig2, ax_arr2 = plt.subplots()
    ax_arr2.scatter(,O32_observed, label= 'O32 Ratio')       
    ax_arr2.set_xlabel(r'Raw log($O_{32}$')
    ax_arr2.set_ylabel(r'Observed log($O_{32}$)')'''


    plt.draw()
    fig.savefig(pdf_pages, format='pdf')
    #print 'end for'
    pdf_pages.close()
    print 'Finished'


def save_for_later():
    lR23 = voronoi[:,2]
    lO32 = voronoi[:,3]
    binnum = voronoi[:,4]
    binnum0= list(set(binnum))
    
    avg_l_R23 = []  #np.zeros(len(binnum0), dtype=np.float64)
    avg_l_O32 = []  #np.zeros(len(binnum0), dtype=np.float64)
    avg_l_O2 = []
    avg_l_O3 = []
    avg_l_Hb = []
    #print type(avg_l_R23), type(avg_l_O32)
    #print avg_l_R23
    for ii in binnum0:
        idx = np.array([xx for xx in range(len(det3)) if binnum[xx] == ii])
            
        avg_O2 = np.mean(O2[idx])/1e-17
        avg_O3 = np.mean(O3[idx])/1e-17
        avg_Hb = np.mean(Hb[idx])/1e-17
        print 'Avg_O3', avg_O3
        list_R23 = np.log10((avg_O2 +avg_O3)/avg_Hb)
        list_O32 = np.log10(avg_O3/avg_O2)
        #print list_R23, type(list_R23), list_O32, type(list_O32)

        avg_l_R23.append(list_R23)
        avg_l_O32.append(list_O32)
        avg_l_O2.append(avg_O2)
        avg_l_O3.append(avg_O3)
        avg_l_Hb.append(avg_Hb)


def comparing_composite_values():
    comparing_composite_values = fitspath +'comparing_composite_values.pdf'
    com_R23 = data1['R_23_Average'].data
    com_O32 = data1['O_32_Average'].data

    OII = data1['OII_3727_Flux_Observed'].data
    OIII4959 = data1['OIII_4958_Flux_Observed'].data
    OIII5007 = data1['OIII_5007_Flux_Observed'].data
    H_BETA = data1['HBETA_Flux_Observed'].data
    binnum = data1['N_Galaxies'].data
    #print binnum

    pdf_pages= PdfPages(comparing_composite_values)
    R23_composite = np.log10((OII + 1.33 *OIII5007)/H_BETA)
    O32_composite = np.log10((1.33*OIII5007)/OII)

    fig, ax_arr = plt.subplots()
    for ii in range(len(binnum)):
        ax_arr.plot(R23_composite[ii], O32_composite[ii], '*')
        print 'calculated', R23_composite[ii], O32_composite[ii]
        ax_arr.plot(com_R23[ii], com_O32[ii], '+')
        print 'computer', com_R23[ii], com_O32[ii]
        ax_arr.set_xlabel('R23')
        ax_arr.set_ylabel('O32')
        fig.savefig(pdf_pages, format='pdf')
    
    pdf_pages.close()
