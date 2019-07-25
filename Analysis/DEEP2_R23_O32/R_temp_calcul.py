
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
from astropy import units as u
from chun_codes.cardelli import *

fitspath_ini='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

#Constants

a = 13205
b = 0.92506
c = 0.98062

def R_calculation(OIII4363, OIII5007, OIII4959, EBV, k_4636, k_5007):
  
    R_value = OIII4363/(OIII5007+OIII4959)* 10**(0.4*EBV*(k_4363-k_5007))
    return R_value  

def temp_calculation(R):
    #T_e = a(-log(R)-b)^(-c)
    T_e =  a*(-np.log10(R)-b)**(-1*c)     
    print T_e  
    return T_e  


def metalicity_calculation(T_e,der_3727_HBETA, der_4959_HBETA, der_5007_HBETA):                 #OIII5007, OIII4959, OIII4363, HBETA, OII3727):
    #12 +log(O+/H) = log(OII/Hb) +5.961 +1.676/t_2 - 0.4logt_2 - 0.034t_2 + log(1+1.35x)
    #12 +log(O++/H) = log(OIII/Hb)+6.200+1.251/t_3 - 0.55log(t_3) - 0.014(t_3)
    #t_2 = 0.7*t_3 +0.17
    
    two_beta = der_3727_HBETA                     #OII3727/HBETA
    three_beta= der_4959_HBETA+ der_5007_HBETA    #(OIII4959+OIII5007)/HBETA
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
    
def run_function(fitspath, dataset, out_ascii, out_fits, pdf_name,  combine_flux_ascii, dustatt= False):  #combine_fits, header
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

    SN_Hgamma     = combine_fits['Hgamma_S/N'].data
    SN_5007       = combine_fits['OIII_5007_S/N'].data
    SN_4959       = combine_fits['OIII_4958_S/N'].data
    SN_4363       = combine_fits['OIII_4363_S/N'].data
    SN_HBETA      = combine_fits['HBETA_S/N'].data
    SN_3727       = combine_fits['OII_3727_S/N'].data

    #Dust Attenuation Values Call
    ###Every time I change the binning method I need to go back and recalculate the dust attentuation##
    print 'Using attenuated dust values'
    dust = '/Users/reagenleimbach/Desktop/Zcalbase_gal/dust_attentuation_values.tbl'
    atten_val = asc.read(dust)

    k_3727 = atten_val['A_3727']
    k_HDELTA = atten_val['A_HDELTA']
    k_Hgamma = atten_val['A_Hgamma']
    k_HBETA = atten_val['A_HBETA']
    k_4363 = atten_val['A_4363']
    k_4959 = atten_val['A_4959']
    k_5007 = atten_val['A_5007']
    
    if dustatt == False: EBV = np.zeros(len(ID))
    if dustatt == True: EBV = atten_val['E(B-V)']


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
    3727_HBETA = OII3727/HBETA
    5007_HBETA = OIII5007/HBETA
    4959_HBETA = OIII4959/HBETA
    4363_5007  = OIII4363/OIII5007
    4363_4959  = OIII4363/OIII4959
    
    5007_3727  = OIII5007/OII3727
    4959_3727  = OIII4959/OII3727

    #Attenuated Ratios
    der_4363_5007  = 4363_5007 * 10**(0.4*EBV*(k_4363-k_5007))
    der_4363_4959  = 4363_4959 * 10**(0.4*EBV*(k_4363-k_4959))
    der_3727_HBETA = 3727_HBETA* 10**(0.4*EBV*(k_3727-k_HBETA))
    der_4959_HBETA = 4959_HBETA* 10**(0.4*EBV*(k_4959-k_HBETA))
    der_5007_HBETA = 5007_HBETA* 10**(0.4*EBV*(k_5007-k_HBETA))
    
    
    #Raw Data
    R_value= R_calculation(OIII4363, OIII5007, OIII4959,EBV, k_4363, k_5007)   #, SN_4636, SN_5007, SN_495)
    T_e= temp_calculation(R_value)  #, R_std)
    O_s_ion, O_d_ion, com_O_log, log_O_s, log_O_d = metalicity_calculation(T_e,der_3727_HBETA, der_4959_HBETA, der_5007_HBETA)  #OIII5007, OIII4959, OIII4363, HBETA, OII3727)


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
    Te_marker = []
    #print Te_marker
    #print('!!!!!', len(OIII4363))
    for oo in range(len(OIII4363)):
        #print indicate[oo]
        #print mDect[oo]
        if indicate[oo]== 0 and mDect[oo]== 0: Marker= '<'
        if indicate[oo]== 1 and mDect[oo]== 0: Marker= '<'
        if indicate[oo]== 0 and mDect[oo]== 1:  Marker = '.'
        if indicate[oo]== 1 and mDect[oo]== 1:  Marker = '.'
        Te_marker.append(Marker)
    print Te_marker

    M_marker = []
    print M_marker
    for ww in range(len(OIII4363)):
        #print indicate[ww]
        #print mDect[ww]
        if indicate[ww]== 0 and mDect[ww]== 0: Marker= '^'
        if indicate[ww]== 1 and mDect[ww]== 0: Marker= '^'
        if indicate[ww]== 0 and mDect[ww]== 1:  Marker = '.'
        if indicate[ww]== 1 and mDect[ww]== 1:  Marker = '.'
        M_marker.append(Marker)
    #print M_marker


    
        
    
    fig1, ax1 = plt.subplots()
    
    for aa in range(len(ID)):
        ax1.scatter(lTe[aa], R23_composite[aa], marker = Te_marker[aa], color = 'b')
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
        ax5.scatter(T_e[a], R23_composite[a], marker = Te_marker[a], color = 'b')
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
        ax2.scatter(lTe[cc], O32_composite[cc], marker = Te_marker[cc], color = 'b')
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
        ax6.scatter(T_e[c], O32_composite[c], marker = Te_marker[c], color = 'b')
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
        ax3.scatter(R23_composite[zz], com_O_log[zz], marker = M_marker[zz], color = 'b')
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
        ax4.scatter(O32_composite[ww], com_O_log[ww], marker = M_marker[ww], color = 'b')
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
        ax7.scatter(x_value[r], y_value[r], marker = Te_marker[r], color = 'b')
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
   


def dust_attenuation(combine_ascii):
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    
    combine_asc = asc.read(combine_ascii)
    ini_con = 0.468
    ID = combine_asc['ID']
    HBeta= combine_asc['HBETA_Flux_Observed'].data
    HGamma= combine_asc['Hgamma_Flux_Observed'].data
    
    lam0_OII = combine_asc['OII_3727_X_bar'].data
    lam0_HDELTA = combine_asc['HDELTA_X_bar'].data
    lam0_Hgamma = combine_asc['Hgamma_X_bar'].data
    lam0_HBETA = combine_asc['HBETA_X_bar'].data
    lam0_4363 = combine_asc['OIII_4363_X_bar'].data
    lam0_4958 = combine_asc['OIII_4958_X_bar'].data
    lam0_5007 = combine_asc['OIII_5007_X_bar'].data

    k_3727 = call_cardelli(lam0_OII)
    k_HDELTA = call_cardelli(lam0_HDELTA)
    k_Hgamma = call_cardelli(lam0_Hgamma)
    k_HBETA = call_cardelli(lam0_HBETA)
    k_4363 = call_cardelli(lam0_4363)
    k_4958 = call_cardelli(lam0_4958)
    k_5007 = call_cardelli(lam0_5007)

    
    
    EBV= np.log10((HBeta/HGamma)*(ini_con))*2.5*(1/(k_Hgamma-k_HBETA))
    for nn in range(len(HGamma)):
        if EBV[nn] <= 0: EBV[nn] = 0
    
    
    A_3727 = EBV*k_3727
    A_HDELTA = EBV*k_HDELTA
    A_Hgamma = EBV*k_Hgamma
    A_HBETA = EBV*k_HBETA
    A_4363 = EBV*k_4363
    A_4958 = EBV*k_4958
    A_5007 = EBV*k_5007


    out_ascii = fitspath_ini+'/dust_attentuation_values.tbl'
    #if not exists(out_ascii_single):
    n2= ('ID','k_3727', 'k_HDELTA', 'k_Hgamma', 'k_HBETA', 'k_4363', 'k_4958', 'k_5007', 'E(B-V)')
    tab1 = Table([ID,k_3727, k_HDELTA, k_Hgamma, k_HBETA , k_4363, k_4958, k_5007, EBV], names=n2)
    asc.write(tab1, out_ascii, format='fixed_width_two_line')
    
    
def call_cardelli(lam0): #, extrapolate=False):
    #lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]* u.angstrom
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    lambda0 = lam0*u.angstrom
    k_values= cardelli(lambda0,R=3.1)
    return k_values


'''k_ascii = fitspath_ini+'/cardelli_k_values.tbl'
    n3 = ('Line_Name','K_value')
    tab3 = Table([line_name, k_values])
    asc.write(tab3, k_ascii, format = 'fixed_width_two_line')'''
