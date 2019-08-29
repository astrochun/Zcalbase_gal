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
import glob

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'


from Zcalbase_gal.Analysis.DEEP2_R23_O32 import general

def all_spectra(): 
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data = general.get_det3()

    x= np.log10(R23)
    y= np.log10(O32)
    fig,ax = plt.subplots()
    ax.scatter(x,SNR3, 1.5)
    ax.scatter(y,SNR3, 1.5)
    plt.show()


###For each Bin first draft###
def SN_each_bin_firstdraft(fitspath, name, asc_tab):
    #fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/'
    #name = 'SN_each_bin.pdf'
    #acs_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Double_Bin_2d_binning_datadet3.tbl'

    pdfpages = PdfPages(fitspath+name)
    ascii_tab = asc.read(asc_tab)
    R23_ini= ascii_tab['R23']
    O32_ini = ascii_tab['O32']
    SN_5007 = ascii_tab['SN_5007']
    N_bin = ascii_tab['N_bin']

    #for ii in range(len(N_bin)):
        #g1 =  np.where((N_bin == ii) & (N_bin ==(ii+1)))[0]

    g1 = np.where((N_bin == 0) |(N_bin ==1))[0]
    g2 = np.where((N_bin == 2) | (N_bin ==3))[0]
    g3 = np.where((N_bin == 4) | (N_bin ==5))[0]
    g4 = np.where((N_bin == 6) | (N_bin ==7))[0]
    g5 = np.where((N_bin == 9) | (N_bin ==8))[0]
    g6 = np.where((N_bin == 11) | (N_bin ==10))[0]
    g7 = np.where((N_bin == 13) | (N_bin ==12))[0]
    g8 = np.where((N_bin == 15) | (N_bin ==14))[0]
    g9 = np.where((N_bin == 17) | (N_bin ==16))[0]

    c1 = np.where((N_bin ==0))[0]
    c2 = np.where((N_bin ==2))[0]
    c3 = np.where((N_bin ==4))[0]
    c4 = np.where((N_bin ==6))[0]
    c5 = np.where((N_bin ==8))[0]
    c6 = np.where((N_bin ==10))[0]
    c7 = np.where((N_bin ==12))[0]
    c8 = np.where((N_bin ==14))[0]
    c9 = np.where((N_bin ==16))[0]
   
    
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    fig6, ax6 = plt.subplots()
    fig7, ax7 = plt.subplots()
    fig8, ax8 = plt.subplots()
    fig9, ax9 = plt.subplots()

    R23 = np.log10(R23_ini)
    O32 = np.log10(O32_ini)
    #SN_5007 = np.log10(SN_5007_ini)
    
    ax1.scatter(R23[g1],SN_5007[g1],1.5, color = 'r', label = 'High Ionization')
    ax1.set_xlabel('Log(R23)')
    ax1.set_ylabel('(SN_5007)')
    ax1.scatter(R23[c1],SN_5007[c1],1.5, color ='b', label = 'Low Ionization')
    ax1.legend(loc=0, fontsize = 3)
    ax1.set_yscale("log")
    fig1.savefig(pdfpages, format='pdf')
    fig1.clear()

    ax2.scatter(R23[g2],SN_5007[g2],1.5, color = 'r', label = 'High Ionization')
    #x = np.median(O32[g2])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax2.set_xlabel('Log(R23)')
    ax2.set_ylabel('(SN_5007)')
    ax2.scatter(R23[c2],SN_5007[c2],1.5, color ='b')
    ax2.legend(loc=0, fontsize = 3)
    ax2.set_yscale("log")
    fig2.savefig(pdfpages, format='pdf')
    fig2.clear()


    ax3.scatter(R23[g3],SN_5007[g3],1.5, color = 'r', label = 'High Ionization')
    #x = np.median(O32[g3])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax3.set_xlabel('Log(R23)')
    ax3.set_ylabel('(SN_5007)')
   
    ax3.scatter(R23[c3],SN_5007[c3],1.5, color ='b', label = 'Low Ionization')
    ax3.set_yscale("log")
    ax3.legend(loc=0, fontsize = 3)
    fig3.savefig(pdfpages, format='pdf')
    fig3.clear()


    ax4.scatter(R23[g4],SN_5007[g4],1.5, color = 'r', label = 'High Ionization')
    #x = np.median(O32[g4])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax4.set_xlabel('Log(R23)')
    ax4.set_ylabel('(SN_5007)')
    
    ax4.scatter(R23[c4],SN_5007[c4],1.5, color = 'b', label = 'Low Ionization')
    ax4.legend(loc=0, fontsize = 3)
    ax4.set_yscale("log")
    fig4.savefig(pdfpages, format='pdf')
    fig4.clear()


    ax5.scatter(R23[g5],SN_5007[g5],1.5, color = 'r', label = 'High Ionization')
    #x = np.median(O32[g5])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax5.set_xlabel('Log(R23)')
    ax5.set_ylabel('(SN_5007)')

    ax5.scatter(R23[c5],SN_5007[c5],1.5, color = 'b', label = 'Low Ionization')
    ax5.legend(loc=0, fontsize = 3)
    ax5.set_yscale("log")
    fig5.savefig(pdfpages, format='pdf')
    fig5.clear()


    ax6.scatter(R23[g6],SN_5007[g6],1.5, color = 'r', label = 'High Ionization')
    #x = np.median(O32[g6])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax6.set_xlabel('Log(R23)')
    ax6.set_ylabel('(SN_5007)')
    
    ax6.scatter(R23[c6],SN_5007[c6],1.5, color ='b', label = 'Low Ionization')
    ax6.legend(loc=0, fontsize = 3)
    ax6.set_yscale("log")
    fig6.savefig(pdfpages, format='pdf')
    fig6.clear()


    ax7.scatter(R23[g7],SN_5007[g7],1.5, color = 'r', label = 'High Ionization')
    
    #x = np.median(O32[g7])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax7.set_xlabel('Log(R23)')
    ax7.set_ylabel('(SN_5007)')
    ax7.scatter(R23[c7],SN_5007[c7],1.5, color = 'b', label = 'Low Ionization')
    ax7.legend(loc=0, fontsize = 3)
    ax7.set_yscale("log")
    fig7.savefig(pdfpages, format='pdf')
    fig7.clear()


    ax8.scatter(R23[g8],SN_5007[g8],1.5, color ='r', label = 'High Ionization')

    #x = np.median(O32[g8])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax8.set_xlabel('Log(R23)')
    ax8.set_ylabel('(SN_5007)')
    ax8.scatter(R23[c8],SN_5007[c8],1.5, color ='b', label = 'Low Ionization')
    ax8.legend(loc=0, fontsize = 3)
    ax8.set_yscale("log")
    fig8.savefig(pdfpages, format='pdf')
    fig8.clear()

    ax9.scatter(R23[g9],SN_5007[g9],1.5, color= 'r', label = 'High Ionization')
    #x = np.median(O32[g9])
    #plt.axvline(x = x, linewidth= 0.3, color= 'k')
    ax9.set_xlabel('Log(R23)')
    ax9.set_ylabel('(SN_5007)')
    
    ax9.scatter(R23[c9],SN_5007[c9],1.5, color= 'b', label = 'Low Ionization')
    ax9.legend(loc=0, fontsize = 3)
    ax9.set_yscale("log")
    fig9.savefig(pdfpages, format='pdf')
    fig9.clear()

    
    pdfpages.close()

def SN_each_bin(fitspath, name, asc_tab):
    #average, median, min, and max S/N
    #fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/'
    #name = 'SN_each_bin.pdf'
    #acs_tab = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_att_200/Double_Bin_2d_binning_datadet3.tbl'


    pdfpages = PdfPages(fitspath+name)
    ascii_tab = asc.read(asc_tab)
    R23_ini= ascii_tab['R23']
    O32_ini = ascii_tab['O32']
    SN_5007 = ascii_tab['SN_5007']
    Bin_number = ascii_tab['Bin_number']
    
    

    R23 = np.log10(R23_ini)
    O32 = np.log10(O32_ini)

    Bin = range(0,27,1)
    all_average = np.zeros(len(Bin))
    all_median = np.zeros(len(Bin))
    all_min = np.zeros(len(Bin))
    all_max = np.zeros(len(Bin))
    
    for ii in range(0,27,1):
        fig1, ax1 = plt.subplots()
        spectra = np.where((Bin_number == ii))[0]
        print ii
        print spectra

        average = np.average(SN_5007[spectra])
        median = np.median(SN_5007[spectra])
        bin_min = np.min(SN_5007[spectra])
        bin_max = np.max(SN_5007[spectra])

        all_average[ii] = average
        all_median[ii] = median
        all_min[ii] = bin_min
        all_max[ii] = bin_max
        
        ax1.scatter(R23[spectra],SN_5007[spectra],1.5, color = 'k', label= "Bin %d"%ii)
        plt.axhline(y= average, linewidth= 0.3, color= 'b', label = 'Average')
        plt.axhline(y= median, linewidth= 0.3, color= 'g', label = 'Median')
        plt.axhline(y= bin_min, linewidth= 0.3, color= 'm', label = 'Minimum')
        plt.axhline(y= bin_max, linewidth= 0.3, color= 'r', label = 'Maximum')

        ax1.set_xlabel('Log(R23)')
        ax1.set_ylabel('SN_5007')
        ax1.legend(loc=0, fontsize = 3)
        ax1.set_yscale("log")
        fig1.savefig(pdfpages, format='pdf')
        fig1.clear()


    n = ('Bin', 'SN_5007 Minimum','SN_5007 Maximum', 'Average', 'Median')
    tab = Table([Bin, all_min, all_max, all_average, all_median], names = n)
    asc.write(tab, fitspath+'/SN_5007_values.tbl', format = 'fixed_width_two_line')
    pdfpages.close()
    
