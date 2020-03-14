#Graphs the temperature, metallicities, R23 and O32 and errors for the individual and composite spectra by importing pre-existing tables and dictionaries


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from scipy.optimize import curve_fit
import scipy.integrate as integ
import glob
from datetime import date

from Metallicity_Stack_Commons import fitspath_reagen as fitspath_ini

def plotting_te_metal(fitspath, out_pdf):
    
    indv_measurements = fitspath + 'individual_derived_properties.tbl'
    composite_file = fitspath + 'bin_derived_properties.tbl'

    #Individual Measurements
    indv_derived = asc.read(indv_measurements)
    comp_derived = asc.read(composite_file)
    
    iID = indv_derived['ID']
    iTe = indv_derived['T_e']
    icom_log = indv_derived['12+log(O/H)']
    ilogR23 = indv_derived['logR23']
    ilogO32 = indv_derived['logO32']
    
    icom_nan = np.isnan(icom_log)
    iidx = np.where((icom_nan ==False))[0]

    iID_idv  = iID[iidx]
    iTe_idv  = iTe[iidx]
    icom_idv = icom_log[iidx]
    iR23_idv = ilogR23[iidx]
    iO32_idv = ilogO32[iidx]
    print("len:", len(iR23_idv))
    
    
    ####DEEP2 and MACT Data#####
    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')

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
    
    #Composite Measurements
    R23_composite=comp_derived['logR23_Composite'].data
    O32_composite=comp_derived['logO32_Composite'].data
    ID_composite=comp_derived['bin_ID'].data
    T_e_composite=comp_derived['T_e'].data
    metal_composite = comp_derived['12+log(O/H)'].data
    ver_detection = comp_derived['Detection'].data

    ver_detect = np.where((ver_detection ==1))[0]
    ver_rlimit = np.where((ver_detection ==0.5))[0]
    nan_detect = np.where((ver_detection == 0))[0]







    
    pdf_pages = PdfPages(fitspath+out_pdf)

    fig, ax = plt.subplots()
    ax.scatter(iR23_idv, iO32_idv, marker = '.', s = 35, color = 'g')
    ax.set_title(r'$R_{23}$ vs. $O_{32}$')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    fig.savefig(pdf_pages, format ='pdf')
    fig.clear()
##################################################################################################
    fig1,ax1 = plt.subplots()
    ax1.scatter(iTe_idv, iR23_idv,  marker = '.', s = 35, color = 'g')
    #ax1.set_title(r'$R_{23}$ vs. $T_e$')
    #ax1.set_xlabel(r'log($R_{23}$)')
    #ax1.set_ylabel('T_e')

    ax1.scatter(T_e_composite[ver_detect], R23_composite[ver_detect], marker = '.',s = 50, color = 'b')
    ax1.scatter(T_e_composite[ver_rlimit], R23_composite[ver_rlimit], marker = '<',s = 35, color = 'b')
    for xx in ver_detect:ax1.annotate(ID_composite[xx], (T_e_composite[xx], R23_composite[xx]), fontsize = '6')
    for xx in ver_rlimit:ax1.annotate(ID_composite[xx], (T_e_composite[xx], R23_composite[xx]), fontsize = '6')
    
    ax1.scatter(der_Te, der_R23, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for b in range(len(ID_der)): ax1.annotate(ID_der[b], (der_Te[b], der_R23[b]), fontsize = '2')

    ax1.scatter(der_Te_MACT, der_R23_MACT, s =20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for q in range(len(ID_der_MACT)): ax1.annotate(ID_der_MACT[q], (der_Te_MACT[q], der_R23_MACT[q]), fontsize= '2')

    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel(r'$R_{23}$')
    ax1.set_title(r'Temperatures vs $R_{23}$ Temperature')




    fig1.savefig(pdf_pages, format ='pdf')
    fig1.clear()


##################################################################################################

    fig2,ax2 = plt.subplots()
    ax2.scatter(iTe_idv,iO32_idv,  marker = '.', s = 35, color = 'g')



    ax2.scatter(T_e_composite[ver_detect], O32_composite[ver_detect], marker = '.',s=50, color = 'b')
    ax2.scatter(T_e_composite[ver_rlimit], O32_composite[ver_rlimit], marker = '<',s=35, color = 'b')
    for c in ver_detect:ax2.annotate(ID_composite[c], (T_e_composite[c], O32_composite[c]), fontsize = '6')
    for c in ver_rlimit:ax2.annotate(ID_composite[c], (T_e_composite[c], O32_composite[c]), fontsize = '6')

    ax2.scatter(der_Te, der_O32, s=20, marker = '*', color = 'k',edgecolors = 'None')
    for f in range(len(ID_der)):ax2.annotate(ID_der[f], (der_Te[f], der_O32[f]), fontsize = '2')

    ax2.scatter(der_Te_MACT, der_O32_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors ='None')
    for s in range(len(ID_der_MACT)):ax2.annotate(ID_der_MACT[s], (der_Te_MACT[s], der_O32_MACT[s]), fontsize= '2')
    
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel(r'$O_{32}$')
    ax2.set_title(r'Temperatures vs $O_{32}$ Temperature')

    fig2.savefig(pdf_pages, format ='pdf')
    fig2.clear()
##################################################################################################
    fig3,ax3 = plt.subplots()
    ax3.scatter(iR23_idv, icom_idv,  marker = '.', s = 35, color = 'g')
    ax3.scatter(R23_composite[ver_detect], metal_composite[ver_detect], marker = '.', s = 50, color = 'b')
    ax3.scatter(R23_composite[ver_rlimit], metal_composite[ver_rlimit], marker = '^', s = 35, color = 'b')
    for zz in ver_detect:ax3.annotate(ID_composite[zz], (R23_composite[zz],metal_composite[zz]), fontsize = '6')
    for zz in ver_rlimit:ax3.annotate(ID_composite[zz], (R23_composite[zz],metal_composite[zz]), fontsize = '6')
    
    ax3.scatter(der_R23, der_OH, s= 20, marker = '*', color = 'k', edgecolors = 'None')
    for gg in range(len(ID_der)): ax3.annotate(ID_der[gg], (der_R23[gg], der_OH[gg]), fontsize = '2')
    ax3.scatter(der_R23_MACT, der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors='None')
    for g in range(len(ID_der_MACT)): ax3.annotate(ID_der_MACT[g], (der_R23_MACT[g], der_OH_MACT[g]), fontsize= '2')

    ax3.set_xlabel(r'$R_{23}$')
    ax3.set_ylabel('12+log(O/H)')
    ax3.set_title(r'$R_{23}$ vs. Composite Metallicity')


    
    fig3.savefig(pdf_pages, format ='pdf')
    fig3.clear()
##################################################################################################
    fig4,ax4 = plt.subplots()
    ax4.scatter(iO32_idv, icom_idv,  marker = '.', s = 35, color = 'g')
    ax4.scatter(O32_composite[ver_detect], metal_composite[ver_detect], marker = '.',s =50, color = 'b')
    ax4.scatter(O32_composite[ver_rlimit], metal_composite[ver_rlimit], marker = '^',s =35, color = 'b')

    for ww in ver_detect:ax4.annotate(ID_composite[ww], (O32_composite[ww], metal_composite[ww]), fontsize = '6')
    for ww in ver_rlimit:ax4.annotate(ID_composite[ww], (O32_composite[ww], metal_composite[ww]), fontsize = '6')
    
    ax4.scatter(der_O32,der_OH, s=20, marker = '*', color = 'k', edgecolors = 'None')
    for hh in range(len(ID_der)): ax4.annotate(ID_der[hh], (der_O32[hh], der_OH[hh]), fontsize = '2')

    ax4.scatter(der_O32_MACT,der_OH_MACT, s=20, marker = 'P', color = 'r', alpha = 0.5, edgecolors = 'None')
    for h in range(len(ID_der_MACT)): ax4.annotate(ID_der_MACT[h], (der_O32_MACT[h], der_OH_MACT[h]), fontsize= '2')

    ax4.set_xlabel(r'$O_{32}$')
    ax4.set_ylabel('12+log(O/H)')
    ax4.set_title(r'$O_{32}$ vs. Composite Metallicity')
    fig4.savefig(pdf_pages, format ='pdf')
    fig4.clear()
        
##################################################################################################        
    pdf_pages.close()
        

