
###THIS CODE IS RUN IN THE GENERAL FUNCTION###

#Calculates the R value, electron temperature, and metallicity from the flux table
#produced by the zoom_and_gauss_general functions

#Currently running: Grid
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
from mpl_toolkits.mplot3d import Axes3D
import sys

#For generalizing for several users
from getpass import getuser

fitspath_ini='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

#Constants

a = 13205
b = 0.92506
c = 0.98062

def R_calculation(OIII4363, OIII5007, OIII4959):
  
    R_value = OIII4363/(OIII5007+OIII4959)
    return R_value  

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)     
    print T_e  
    return T_e  


def metalicity_calculation(T_e,OIII5007, OIII4959, OIII4363, HBETA, OII3727):
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17
    
    two_beta = OII3727/HBETA
    three_beta= (OIII4959+OIII5007)/HBETA
    t_3 = T_e*1e-4
    t_2 = 0.7*t_3 +0.17
    x2 = 1e-4 * 1e3 * t_2**(-0.5)

    O_s_ion_log = np.log10(two_beta) +5.961 +1.676/t_2 - 0.4*np.log10(t_2) - 0.034*t_2 + np.log10(1+1.35*x2)-12
    O_d_ion_log = np.log10(three_beta)+6.200+1.251/t_3 - 0.55*np.log10(t_3) - 0.014*(t_3)-12

    O_s_ion = 10**(O_s_ion_log)
    O_d_ion = 10**(O_d_ion_log)
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) +12

    return O_s_ion , O_d_ion, com_O_log, O_s_ion_log, O_d_ion_log

def limit_function(combine_flux_ascii):
    #hgamma/hgamma_sn*3 = 3sigma
    combine_fits= asc.read(combine_flux_ascii)

    Hgamma = combine_fits['Hgamma_Flux_Observed'].data
    Hgamma_SN    = combine_fits['Hgamma_S/N'].data

    up_temp = (Hgamma/Hgamma_SN) *3

    print 'up_temp', up_temp
    return up_temp
    
def run_function(fitspath, dataset, out_ascii, out_fits, pdf_name,  combine_flux_ascii):  #combine_fits, header
    verification_table = fitspath_ini+'/Master_Verification_Tables/'+dataset+'_table.tbl'
    mver_tab = asc.read(verification_table)


    #Fits Table Calls
    #combine_fits, header = fits.getdata(combine_flux_table, header = True)
    combine_fits= asc.read(combine_flux_ascii)

    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')

    #Ascii Table from FITTING 
    OIII5007 = combine_fits['OIII_5007_Flux_Observed'].data
    OIII4959 = combine_fits['OIII_4958_Flux_Observed'].data
    raw_OIII4363 = combine_fits['OIII_4363_Flux_Observed'].data
    Hgamma = combine_fits['Hgamma_Flux_Observed'].data
    HBETA    = combine_fits['HBETA_Flux_Observed'].data
    OII3727  = combine_fits['OII_3727_Flux_Observed'].data
    R23_avg      = combine_fits['R_23_Average'].data
    O32_avg      = combine_fits['O_32_Average'].data
    N_Galaxy = combine_fits['N_Galaxies'].data
    ID = combine_fits['ID'].data

    SN_Hgamma    = combine_fits['Hgamma_S/N'].data
    SN_5007       = combine_fits['OIII_5007_S/N'].data
    SN_4959       = combine_fits['OIII_4958_S/N'].data
    SN_4363       = combine_fits['OIII_4363_S/N'].data
    SN_HBETA      = combine_fits['HBETA_S/N'].data
    SN_3727       = combine_fits['OII_3727_S/N'].data

    #DEEP2 Derived 
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_Te = derived['Te'].data
    der_OH = derived['OH'].data
    ID_der = derived['ID'].data
    #der_OH_log = np.log10(er_OH_log)

    #MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_Te_MACT = derived_MACT['Te'].data
    der_OH_MACT = derived_MACT['OH'].data
    ID_der_MACT = derived_MACT['ID'].data
    print len(der_Te_MACT), len(der_R23_MACT)

    R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
    O32_composite = np.log10((1.33*OIII5007)/OII3727)

    print 'R23_composite', R23_composite
    print 'O32_composite', O32_composite 
    
    up_limit = (Hgamma/SN_Hgamma) *3
    print 'up_limit', up_limit

    OIII4363 = np.zeros(len(raw_OIII4363))
    indicate = np.zeros(len(raw_OIII4363))
    for ii in range(len(OIII4363)):
        print SN_4363[ii]
        if SN_4363[ii] >= 3:
            print 'regular'
            print '4363:' , raw_OIII4363[ii]
            OIII4363[ii]= raw_OIII4363[ii]
            indicate[ii]= 1 
        else:
            print 'upper limit'
            print '4363: ', up_limit[ii]
            OIII4363[ii]= up_limit[ii]
            indicate[ii]= 0
    print OIII4363
    print indicate 

    #Line Ratios
    '''3727_HBETA = OII3727/HBETA
    5007_HBETA = OIII5007/HBETA
    4959_HBETA = OIII4959/HBETA
    5007_3727  = OIII5007/OII3727
    4959_3727  = OIII4959/OII3727
    4363_5007  = OIII4363/OIII5007'''
    
    #Raw Data
    R_value= R_calculation(OIII4363, OIII5007, OIII4959)   #, SN_4636, SN_5007, SN_495)
    T_e= temp_calculation(R_value)  #, R_std)
    O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metalicity_calculation(T_e,OIII5007, OIII4959, OIII4363, HBETA, OII3727)


    #if not exists(out_ascii):
    n=  ('ID', 'Detection', 'R23_Composite', 'O32_Composite', 'R_23_Average', 'O_32_Average', 'N_Galaxies', 'Observed_Flux_5007', 'S/N_5007', 'Observed_Flux_4959', 'S/N_4959', 'Observed_Flux_4363', 'S/N_4363', 'Observed_Flux_HBETA', 'S/N_HBETA', 'Observed_Flux_3727', 'S/N_3727', 'Temperature', 'O_s_ion', 'O_d_ion', 'com_O_log' )  #, '3727_HBETA', '5007_HBETA', '4959_HBETA', '5007_3727', '4959_3727', '4363_5007')
        
    tab0 = Table([ID , indicate, R23_composite, O32_composite, R23_avg, O32_avg, N_Galaxy, OIII5007, SN_5007, OIII4959, SN_4959, OIII4363, SN_4363, HBETA, SN_HBETA, OII3727, SN_3727, T_e, O_s_ion, O_d_ion, com_O_log], names=n)   # 3727_HBETA,5007_HBETA,4959_HBETA,5007_3727, 4959_3727, 4363_5007 
    asc.write(tab0, out_ascii, format='fixed_width_two_line')

    
    tab0.write(out_fits,format = 'fits', overwrite = True)

    #Plots
    #name = 'Grid_temperature_vs_R23.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)
    color_len = len(R23_composite)
    
    color_arr = plt.cm.get_cmap('rainbow')
    
    print 'T_e', T_e
    print 'R23', R23_composite 
    print 'O32', O32_composite

    lTe = np.log10(T_e)
    lder_Te = np.log10(der_Te)
    lMACT_Te = np.log10(der_Te_MACT)

    mDect = mver_tab['Detection'].data
    marker = []
    print marker
    for oo in range(len(OIII4363)):
        print indicate[oo]
        print mDect[oo]
        if indicate[oo]== 0 and mDect[oo]== 0: Marker= '^'
        if indicate[oo]== 1 and mDect[oo]== 0: Marker= '^'
        if indicate[oo]== 0 and mDect[oo]== 1:  Marker = '.'
        if indicate[oo]== 1 and mDect[oo]== 1:  Marker = '.'
        marker.append( Marker)
    print marker
        
    
    fig1, ax1 = plt.subplots()
    
    for aa in range(len(ID)):
        ax1.scatter(lTe[aa], R23_composite[aa], marker = marker[aa], color = 'b')
        ax1.annotate(ID[aa], (lTe[aa], R23_composite[aa]), fontsize = '6')
    ax1.scatter(lder_Te, der_R23, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for bb in range(len(ID_der)):
        ax1.annotate(ID_der[bb], (lder_Te[bb], der_R23[bb]), fontsize ='2')
    ax1.scatter(lMACT_Te, der_R23_MACT, s = 20, marker = 'P', color = 'r', edgecolors = 'None', alpha = 0.5)
    for qq in range(len(ID_der_MACT)):
        ax1.annotate(ID_der_MACT[qq], (lMACT_Te[qq], der_R23_MACT[qq]), fontsize= '2')
    ax1.set_xlabel('Log Temperature (K)')
    ax1.set_ylabel('R_23')
    ax1.set_title('Temperatures_vs_R23')
    pdf_pages.savefig()

    #With Limits
    fig5, ax5 = plt.subplots()
    
    for a in range(len(ID)):
        ax5.scatter(T_e[a], R23_composite[a], marker = marker[a], color = 'b')
        ax5.annotate(ID[a], (T_e[a], R23_composite[a]), fontsize = '6')
    ax5.scatter(der_Te, der_R23, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for b in range(len(ID_der)):
        ax5.annotate(ID_der[b], (der_Te[b], der_R23[b]), fontsize = '2')
    ax5.scatter(der_Te_MACT, der_R23_MACT, s =20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for q in range(len(ID_der_MACT)):
        ax5.annotate(ID_der_MACT[q], (der_Te_MACT[q], der_R23_MACT[q]), fontsize= '2')
    ax5.set_xlabel('Temperature (K)')
    ax5.set_ylabel('R_23')
    ax5.set_title('Temperatures_vs_R23 with Limits on Temperature')
    ax5.set_xlim(5000,21500)
    pdf_pages.savefig()
     
    fig2, ax2 = plt.subplots()
    
    for cc in range(len(ID)):
        ax2.scatter(lTe[cc], O32_composite[cc], marker = marker[cc], color = 'b')
        ax2.annotate(ID[cc], (lTe[cc], O32_composite[cc]), fontsize = '6')
    ax2.scatter(lder_Te, der_O32, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for ff in range(len(ID_der)):
        ax2.annotate(ID_der[ff], (lder_Te[ff], der_O32[ff]), fontsize = '2')

    ax2.scatter(lMACT_Te, der_O32_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5,edgecolors = 'None')
    for ss in range(len(ID_der_MACT)):
        ax2.annotate(ID_der_MACT[ss], (lMACT_Te[ss], der_O32_MACT[ss]), fontsize= '2')
    ax2.set_xlabel('Log Temperature (K)')
    ax2.set_ylabel('O_32')
    ax2.set_title('Temperatures_vs_O32')
    #ax2.set_xlim(1000,21500)
    pdf_pages.savefig()

    #With Limits
    fig6, ax6 = plt.subplots()
    
    for c in range(len(ID)):
        ax6.scatter(T_e[c], O32_composite[c], marker = marker[c], color = 'b')
        ax6.annotate(ID[c], (T_e[c], O32_composite[c]), fontsize = '6')
    ax6.scatter(der_Te, der_O32, s=20, marker = '*', color = 'k',edgecolors = 'None')
    for f in range(len(ID_der)):
        ax6.annotate(ID_der[f], (der_Te[f], der_O32[f]), fontsize = '2')

    ax6.scatter(der_Te_MACT, der_O32_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors ='None')
    for s in range(len(ID_der_MACT)):
        ax6.annotate(ID_der_MACT[s], (der_Te_MACT[s], der_O32_MACT[s]), fontsize= '2')
    ax6.set_xlabel('Temperature (K)')
    ax6.set_ylabel('O_32')
    ax6.set_title('Temperatures_vs_O32 with Limits on Temperature')
    ax6.set_xlim(5000,21500)
    pdf_pages.savefig()

    fig3, ax3 = plt.subplots()
    
    for zz in range(len(ID)):
        ax3.scatter(R23_composite[zz], com_O_log[zz], marker = marker[zz], color = 'b')
        ax3.annotate(ID[zz], (R23_composite[zz],com_O_log[zz]), fontsize = '6')
    ax3.scatter(der_R23, der_OH, s= 20, marker = '*', color = 'k', edgecolors = 'None')
    for gg in range(len(ID_der)):
        ax3.annotate(ID_der[gg], (der_R23[gg], der_OH[gg]), fontsize = '2')
    ax3.scatter(der_R23_MACT, der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors='None')
    for g in range(len(ID_der_MACT)):
        ax3.annotate(ID_der_MACT[g], (der_R23_MACT[g], der_OH_MACT[g]), fontsize= '2')
    ax3.set_xlabel('R23')
    ax3.set_ylabel('12+log(O/H) Te')
    ax3.set_title('R23 vs. Composite Metalicity')
    #ax3.plot(BR23,B_com_R23, 'k')
    #ax2.set_xlim(1000,21500)
    pdf_pages.savefig()

    fig4, ax4 = plt.subplots()
    
    for ww in range(len(ID)):
        ax4.scatter(O32_composite[ww], com_O_log[ww], marker = marker[ww], color = 'b')
        ax4.annotate(ID[ww], (O32_composite[ww], com_O_log[ww]), fontsize = '6')
    ax4.scatter(der_O32,der_OH, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for hh in range(len(ID_der)):
        ax4.annotate(ID_der[hh], (der_O32[hh], der_OH[hh]), fontsize = '2')
    ax4.scatter(der_O32_MACT,der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for h in range(len(ID_der_MACT)):
        ax4.annotate(ID_der_MACT[h], (der_O32_MACT[h], der_OH_MACT[h]), fontsize= '2')
    ax4.set_xlabel('O32')
    ax4.set_ylabel('12+log(O/H) Te')
    ax4.set_title('O32 vs. Composite Metalicity')
    #ax4.plot(BO32, B_com_O32, 'k')
    #ax2.set_xlim(1000,21500)
    pdf_pages.savefig()


    #Can you plot OIII4363/OIII5007 vs T_e and overlay the Te--4363/5007 line ratio solution?  Since T_e is so wide, do a logarithmic plot on that axis.

    fig7, ax7 = plt.subplots()
    x_value = np.log10(OIII4363/OIII5007)
    y_value = np.log10(T_e)
    z_value = np.log10(der_Te)
    
    for r in range(len(ID)):
        ax7.scatter(x_value[r], y_value[r], marker = marker[r], color = 'b')
        ax7.annotate(ID[r], (x_value[r], y_value[r]))
    xxx = np.arange(min(x_value),max(x_value), 1)
    xx_value = np.zeros(len(xxx))
    
    t_E_value = np.zeros(len(xxx))
    for xx in range(len(xxx)):
        xx_value[xx] = xxx[xx]
        T_e =  a*(-np.log10(xxx[xx])-b)**(-1*c)
        t_E_value[xx] = T_e
    ax7.plot(xx_value, t_E_value, 'k')
    
    ax7.set_xlabel('log(OIII4363/OIII5007)')
    ax7.set_ylabel('log(T_e)')
    ax7.set_title('log(OIII4363/OIII5007) vs log(T_e)')
    pdf_pages.savefig()

    pdf_pages.close()

    

    #3D plots
    fig_3d= plt.figure(figsize=(10,8))
    ax_3d = plt.axes(projection='3d')
    ax_3d.set_xlabel('R_23')
    ax_3d.set_ylabel('O_32')
    ax_3d.set_zlabel('Temperature')
    ax_3d.set_zlim(4000,26000)
    ax_3d.scatter(R23_composite, O32_composite, T_e, marker='.', linewidths = None)
    #plt.show()





    fig1.clear()
    fig2.clear()
    fig3.clear()
    fig4.clear()
    fig5.clear()
    fig6.clear()
    fig7.clear()
    fig_3d.clear()
   
