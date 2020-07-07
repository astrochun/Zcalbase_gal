#Graphs the temperature, metallicities, R23 and O32 and errors for the individual and composite spectra by importing pre-existing tables and dictionaries

"""
Keywords:
        fitspath -> path to where files come and are saved to
         revised  -> refers to if using the bin_derived_prop_revised temperature
                     and metallicity measurements which right now implement dust attenuation
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import join
from scipy.optimize import curve_fit
import scipy.integrate as integ
import glob

from Metallicity_Stack_Commons.Metallicity_Stack_Commons import fitspath_reagen as fitspath_ini
from Metallicity_Stack_Commons.Metallicity_Stack_Commons.column_names import filename_dict

def plotting_te_metal(fitspath, revised=False):
    
    indv_all_file = join(fitspath, filename_dict['indv_bin_info'])

    if revised:
        composite_file = join(fitspath, filename_dict['bin_derived_prop_rev'])
        out_pdf = join(fitspath, 'temperature_metallicity_revised.pdf')
    else:
        composite_file = join(fitspath, filename_dict['bin_derived_prop'])
        out_pdf = join(fitspath, 'temperature_metallicity_plots.pdf')

    
    # Individual Measurements
    comp_derived = asc.read(composite_file)
    indv_all = asc.read(indv_all_file)
    
    iID = indv_derived['ID']
    iTe = indv_derived['T_e']
    icom_log = indv_derived['12+log(O/H)']
    ilogR23 = indv_all['logR23']
    ilogO32 = indv_all['logO32']
    print(ilogR23)
    
    #icom_nan = np.isnan(icom_log)
    iidx = np.where((icom_log != 0.0))[0]

    iID_idv  = iID[iidx]
    iTe_idv  = iTe[iidx]
    icom_idv = icom_log[iidx]
    iiR23_idv = ilogR23[iidx]
    iiO32_idv = ilogO32[iidx]
    iR23_idv = np.log10(iiR23_idv)
    iO32_idv = np.log10(iiO32_idv)
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


    pdf_pages = PdfPages(out_pdf)

    fig, ax = plt.subplots()
    ax.scatter(iR23_idv, iO32_idv, marker='.', s=35, color='g')
    ax.set_title(r'$R_{23}$ vs. $O_{32}$')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    fig.savefig(pdf_pages, format ='pdf')
    fig.clear()
##################################################################################################
    fig1, ax1 = plt.subplots()

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
    ax3.scatter(R23_composite[ver_detect], metal_composite[ver_detect], marker = '.', s = 50, color = 'b')
    ax3.scatter(R23_composite[ver_rlimit], metal_composite[ver_rlimit], marker = '^', s = 35, color = 'b')
    for zz in ver_detect:ax3.annotate(ID_composite[zz], (R23_composite[zz],metal_composite[zz]), fontsize = '6')
    for zz in ver_rlimit:ax3.annotate(ID_composite[zz], (R23_composite[zz],metal_composite[zz]), fontsize = '6')
    
    ax3.scatter(der_R23, der_OH, s= 20, marker = '*', color = 'k', edgecolors = 'None')
    for gg in range(len(ID_der)):
        ax3.annotate(ID_der[gg], (der_R23[gg], der_OH[gg]), fontsize='2')
    ax3.scatter(der_R23_MACT, der_OH_MACT, s=20, marker='P', color='r', alpha=0.5, edgecolors='None')
    for g in range(len(ID_der_MACT)):
        ax3.annotate(ID_der_MACT[g], (der_R23_MACT[g], der_OH_MACT[g]), fontsize='2')
    ax3.set_xlim(0.5, 1.1)
    ax3.set_ylim(6.75, 9.25)
    ax3.set_xlabel(r'$R_{23}$')
    ax3.set_ylabel('12+log(O/H)')
    ax3.set_title(r'$R_{23}$ vs. Composite Metallicity')

    fig3.savefig(pdf_pages, format ='pdf')
    fig3.clear()
##################################################################################################
    fig4,ax4 = plt.subplots()
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
        

def Jiang_comparison():
    #log(R23) = a +bx+cx^2 - d(e+x)y
    #x = 12+log(O/H)
    #y = log(O32)
    
    fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/'
    validation = asc.read(join(fitspath, 'bin_validation.revised.tbl'))
    temp_tab = asc.read(join(fitspath, 'bin_derived_properties.tbl'))


    pdf_pages = PdfPages(join(fitspath, 'comparsion_Jiang_Zcal.pdf'))

    bin_ID = temp_tab['bin_ID'].data
    lR23_comp = temp_tab['logR23_Composite'].data
    lO32_comp = temp_tab['logO32_Composite'].data
    zmetal = temp_tab['12+log(O/H)'].data

    valid = validation['Detection'].data
    detect = np.where(valid == 1.0)[0]
    rlimit = np.where(valid == 0.5)[0]

    valid_ID = bin_ID[detect]

    lR23 = lR23_comp[detect]
    lO32 = lO32_comp[detect]

    rlR23 = lR23_comp[rlimit]
    rlO32 = lO32_comp[rlimit]

    metal_det = zmetal[detect]
    metal_rl = zmetal[rlimit]

    derived = asc.read(join(fitspath_ini, 'DEEP2_R23_O32_derived.tbl'))
    derived_MACT = asc.read(join(fitspath_ini, 'MACT_R23_O32_derived.tbl'))
    
    #DEEP2 Derived 
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_OH = derived['OH'].data
    
    #MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_MACT['OH'].data
    
    
    a = -24.135
    b = 6.1532
    c = -0.37866
    d = -0.147
    e = -7.071

    jR23_det = np.zeros(len(metal_det))
    A_comparison = np.zeros(len(metal_det))
    count = len(metal_det) + len(der_R23) + len(der_R23_MACT)
    print(count)

    for ii in range(len(metal_det)):
        jR23_det[ii] = a +b*metal_det[ii]+c*(metal_det[ii]*metal_det[ii]) - d*(e+metal_det[ii])*lO32[ii]
        A_comparison[ii]= (jR23_det[ii] - lR23[ii])
    
    jR23_DEEP = np.zeros(len(der_R23))
    jR23_MACT = np.zeros(len(der_R23_MACT))
    B_comparison = np.zeros(len(der_R23))
    C_comparison = np.zeros(len(der_R23_MACT))
    for bb in range(len(jR23_DEEP)):
        jR23_DEEP[bb] = a +b*der_OH[bb]+c*(der_OH[bb]*der_OH[bb]) - d*(e+der_OH[bb])*der_O32[bb]
        B_comparison[bb]= (jR23_DEEP[bb] - der_R23[bb])
        
    
    for aa in range(len(jR23_MACT)):
        jR23_MACT[aa] = a +b*der_OH_MACT[aa]+c*(der_OH_MACT[aa]*der_OH_MACT[aa]) - d*(e+der_OH_MACT[aa])*der_O32_MACT[aa]
        C_comparison[aa]= (jR23_MACT[aa] - der_R23_MACT[aa])

    arr_sum = np.concatenate((A_comparison, B_comparison, C_comparison), axis=None)
    print(arr_sum)
    med0 = np.median(arr_sum)
    avg0 = np.average(arr_sum)
    sig0 = np.std(arr_sum)
    print('med: ', med0, 'avg: ', avg0, 'sig: ', sig0)
    
    fig, ax = plt.subplots()
    ax.scatter(lR23, jR23_det, marker='D', color='b', alpha=0.75, label='Composite Detections')
    for aa in range(len(valid_ID)):
        ax.annotate(valid_ID[aa], (lR23[aa], jR23_det[aa]), fontsize='6')
    ax.scatter(der_R23, jR23_DEEP, marker='3', color='r', label='DEEP2 Individual Spectra')
    ax.scatter(der_R23_MACT, jR23_MACT, marker='4', color='m', label='MACT Individual Spectra')
    ax.set_xlabel(r'Observed $log(R_{23})$')
    ax.set_ylabel(r'Estimated $log(R_{23})$')
    plt.plot(lR23, lR23, 'k', label='One to one line')
    plt.legend()

    an_txt  = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    #ax.annotate(an_txt, [0.2,0.015], xycoords='axes fraction', va='bottom', ha='right',fontsize=10)
    
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()


def dm_Jiang_comparison():
    fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'
    fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/'
    pdf_pages = PdfPages(fitspath+'DEEPMACT_Jiang_Zcal.pdf')
    
    
    a = -24.135
    b = 6.1532
    c = -0.37866
    d = -0.147
    e = -7.071

    

    fig, ax = plt.subplots()
    
    #for aa in range(len(valid_ID)): ax.annotate(valid_ID[aa], (lR23[aa], jR23_det[aa]), fontsize = '6')
    ax.set_xlabel(r'$R_{23}$ Zcalbase')
    ax.set_ylabel(r'$R_{23}$ Jiang')
    plt.legend()
    plt.plot(der_R23_MACT, der_R23_MACT, 'm', label = 'One to one line')
    
    fig.savefig(pdf_pages, format ='pdf')

    pdf_pages.close()

    
def Bian_comparison(fitspath):
    """

    Returns
    -------

    log(R23) = a +bx+cx^2 - d(e+x)y
    x = 12+log(O/H)
    y = log(O32)
    """
    validation = asc.read(fitspath + 'bin_validation.revised.tbl')
    temp_tab = asc.read(fitspath + 'bin_derived_properties.tbl')


    pdf_pages = PdfPages(fitspath+'comparsion_Bian_Zcal.pdf')

    
    bin_ID = temp_tab['bin_ID']
    lR23_comp = temp_tab['logR23_Composite']
    lO32_comp = temp_tab['logO32_Composite']
    zmetal = temp_tab['12+log(O/H)']

    valid = validation['Detection']
    detect = np.where((valid ==1.0))[0]
    rlimit = np.where((valid ==0.5))[0]

    valid_ID = bin_ID[detect]

    lR23 = lR23_comp[detect]
    lO32 = lO32_comp[detect]

    
    metal_det = zmetal[detect]

    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')
    
    #DEEP2 Derived 
    er_R23 = derived['R23'].data
    er_O32 = derived['O32'].data
    der_R23 = np.log10(er_R23)
    der_O32 = np.log10(er_O32)
    der_OH = derived['OH'].data
    
    #MACT Derived
    er_R23_MACT = derived_MACT['R23'].data
    er_O32_MACT = derived_MACT['O32'].data
    der_R23_MACT = np.log10(er_R23_MACT)
    der_O32_MACT = np.log10(er_O32_MACT)
    der_OH_MACT = derived_MACT['OH'].data

    bR23_DEEP = np.zeros(len(der_R23))
    B_comparison = np.zeros(len(der_R23))
    BO_comparison = np.zeros(len(der_R23))
    bR23_MACT = np.zeros(len(der_R23_MACT))
    C_comparison = np.zeros(len(der_R23_MACT))
    CO_comparison = np.zeros(len(der_R23_MACT))

    bO32_DEEP = np.zeros(len(der_R23))
    bO32_MACT = np.zeros(len(der_R23_MACT))
    

    jR23_det = np.zeros(len(metal_det))
    jO32_det = np.zeros(len(metal_det))
    A_comparison = np.zeros(len(metal_det))
    AO_comparison = np.zeros(len(metal_det))
    for ii in range(len(metal_det)):
        jR23_det[ii] = 138.0430 - 54.8284*metal_det[ii] + 7.2954*metal_det[ii]*metal_det[ii] - 0.32293*metal_det[ii]*metal_det[ii]*metal_det[ii]
        A_comparison[ii] = (jR23_det[ii] - lR23[ii])
        jO32_det[ii] = (-1/0.59)*(metal_det[ii]-8.54)
        AO_comparison[ii] =  (jO32_det[ii] - lO32[ii])
    #for aa in range(len(metal_rl)):
        #jR23_rl[aa] = a +b*metal_rl[aa]+c*(metal_rl[aa]*metal_rl[aa]) - d*(e+metal_rl[aa])*+rlO32[aa]


    for bb in range(len(bR23_DEEP)):
        bR23_DEEP[bb] = 138.0430 - 54.8284* der_OH[bb] + 7.2954* der_OH[bb]* der_OH[bb] - 0.32293* der_OH[bb]* der_OH[bb]* der_OH[bb]
        B_comparison[bb] = (bR23_DEEP[bb] - der_R23[bb])
        bO32_DEEP[bb] = (-1/0.59)*(der_OH[bb]-8.54)
        BO_comparison[bb] = (bO32_DEEP[bb] - der_O32[bb])
    for aa in range(len(bR23_MACT)):
       
        bR23_MACT[aa] = 138.0430 - 54.8284* der_OH_MACT[aa] + 7.2954* der_OH_MACT[aa]* der_OH_MACT[aa] - 0.32293* der_OH_MACT[aa]* der_OH_MACT[aa]* der_OH_MACT[aa]
        C_comparison[aa] = (bR23_MACT[aa] - der_R23_MACT[aa])
        bO32_MACT[aa] = (-1/0.59)*(der_OH_MACT[aa]-8.54)
        CO_comparison[aa] = (bO32_MACT[aa] - der_O32_MACT[aa])


    arr_sum = np.concatenate((A_comparison, B_comparison, C_comparison), axis= None)
    med0 = np.median(arr_sum)
    avg0 = np.average(arr_sum)
    sig0 = np.std(arr_sum)

    print('DEEP2 x: ', der_R23)
    print('DEEP2 y: ', bR23_DEEP)
    print('MACT x:', der_R23_MACT)
    print('MACT y:' ,bR23_MACT)
    print('concatenate array: ', arr_sum)
    print('med: ', med0, 'avg: ', avg0, 'sig: ', sig0)

    n = ('DEEP2 x','DEEP2 y')
    np.savez('xandy.npz', DEEPx = der_R23, DEEPy = bR23_DEEP, MACTx =  der_R23_MACT, MACTy=bR23_MACT)
    n2 = ('MACT x','MACT y')

    fig, ax = plt.subplots()
    ax.scatter(lR23, jR23_det, marker = 'D', color = 'b', alpha = 0.75, label = 'Detections')
    ax.scatter(der_R23,bR23_DEEP,  marker = '3', color = 'r', label = 'DEEP2 Individual Spectra')
    
    ax.scatter(der_R23_MACT,bR23_MACT, marker = '4', color = 'm', label = 'MACT Individual Spectra')
    for aa in range(len(valid_ID)): ax.annotate(valid_ID[aa], (lR23[aa], jR23_det[aa]), fontsize = '6')
    ax.set_xlabel('Observed '+r'$log(R_{23})$')
    ax.set_ylabel('Estimated '+r'$log(R_{23})$')
    plt.plot(lR23, lR23, 'k', label = 'One to one line')
    plt.legend()

    an_txt  = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.annotate(an_txt, [0.2,0.85], xycoords='axes fraction', va='bottom', ha='right',
                fontsize=10)
    
    fig.savefig(pdf_pages, format ='pdf')
    fig.clear()

    arr_sum1 = np.concatenate((AO_comparison, BO_comparison, CO_comparison), axis= None)
    medO0 = np.median(arr_sum1)
    avgO0 = np.average(arr_sum1)
    sigO0 = np.std(arr_sum1)

    print(arr_sum1)
    print('med: ', medO0, 'avg: ', avgO0, 'sig: ', sigO0)
    
    fig, ax = plt.subplots()
    ax.scatter(lO32, jO32_det, marker = 'D', color = 'b', label = 'Detections')
    ax.scatter(der_O32,bO32_DEEP, marker = '3', color = 'r', label = 'DEEP2 Individual Spectra')
    ax.scatter(der_O32_MACT,bO32_MACT, marker = '4', color = 'm', label = 'MACT Individual Spectra')
    for aa in range(len(valid_ID)): ax.annotate(valid_ID[aa], (lO32[aa], jO32_det[aa]), fontsize = '6')
    ax.set_xlabel('Observed '+r'$log(O_{32})$ ')
    ax.set_ylabel('Estimated '+r'$log(O_{32})$')
    plt.plot(jO32_det, jO32_det, 'k', label = 'One to one line')
    plt.legend()

    an_txt  = r'$<\Delta_{R_{23}}>$ : %0.2f' % avgO0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % medO0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sigO0
    ax.annotate(an_txt, [0.2,0.85], xycoords='axes fraction', va='bottom', ha='right',
                fontsize=10)
    
    fig.savefig(pdf_pages, format ='pdf')
    
    
    pdf_pages.close()


def dm_Bian_comparison():
    fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'
    fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0417/'
    pdf_pages = PdfPages(fitspath+'DEEPMACT_bian_Zcal.pdf')
    
    fig, ax = plt.subplots()
    ax.scatter(der_R23, bR23_DEEP, marker='*', color='b', label='DEEP2 Individual Spectra')
    ax.scatter(der_R23_MACT, bR23_MACT, marker='*', color='r', label='MACT Individual Spectra')
    plt.plot(der_R23_MACT, der_R23_MACT, 'm', label='One to one line')
    ax.set_ylim(0.65, 1.1)
    ax.set_xlabel(r'$R_{23}$ Zcalbase')
    ax.set_ylabel(r'$R_{23}$ Bian')
    plt.legend()
    
    
    fig.savefig(pdf_pages, format ='pdf')
    fig.clear()

    fig, ax = plt.subplots()
    ax.scatter(der_O32, bO32_DEEP, marker='*', color='b', label='DEEP2 Individual Spectra')
    ax.scatter(der_O32_MACT, bO32_MACT, marker='*', color='r', label='MACT Individual Spectra')
    plt.plot(der_O32_MACT, der_O32_MACT, 'm', label='One to one line')
    ax.set_xlabel(r'$O_{32}$ Zcalbase')
    ax.set_ylabel(r'$O_{32}$ Bian')
    plt.legend()
    
    
    fig.savefig(pdf_pages, format ='pdf')

    pdf_pages.close()

    '''
    up_pred_det = np.zeros(len(lR23))
    low_pred_det = np.zeros(len(lR23))
    up_pred_rl = np.zeros(len(rlR23))
    low_pred_rl = np.zeros(len(rlR23))
    for ii in range(len(detect)):
        A = lR23[ii]-a+d*e*lO32[ii]
        B = b-d*lO32[ii]
        up_pred_det[ii] = (-B +np.sqrt(B*B - 4*c*A))/(2*c)
        low_pred_det[ii] = (-B -np.sqrt(B*B - 4*c*A))/(2*c)

    for aa in range(len(detect)):
        A = rlR23[aa]-a+d*e*rlO32[aa]
        B = b-d*rlO32[aa]
        up_pred_rl[aa] = (-B +np.sqrt(B*B - 4*c*A))/(2*c)
        low_pred_rl[aa] = (-B -np.sqrt(B*B - 4*c*A))/(2*c)'''


    #print(up_pred_det, low_pred_det,  up_pred_rl, low_pred_rl)
   

    '''
    fig, ax = plt.subplots()
    ax.scatter(lR23, up_pred_det, marker = '*', color = 'b', label = 'Upper')
    ax.scatter(lR23, low_pred_det, marker = '*', color = 'g', label = 'Lower')
    ax.scatter(lR23, metal_det, marker = 'o', color = 'r', label = 'Composite')
    fig.savefig(pdf_pages, format ='pdf')
    
    
    ax.scatter(rlR23, up_pred_rl, marker = '*', color = 'b', label = 'Upper')
    ax.scatter(rlR23, low_pred_rl, marker = '*', color = 'g', label = 'Lower')
    ax.scatter(rlR23, metal_rl, marker = 'o', color = 'r', label = 'Composite')
    fig.savefig(pdf_pages, format ='pdf')'''
