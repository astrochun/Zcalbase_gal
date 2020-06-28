# This function runs the entire process start to finish
# This weekend combine the grid and voronoi if statements, voronoi20, log plots for met
# EW values:equival width
### 

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



##########Dictionary of pdf and other files specific to the Zcalbase_gal project###########
name_dict = dict()

name_dict['gridnpz']  = 'grid.npz'
name_dict['gridpdf']  = 'grid.pdf'
name_dict['Stackname']= 'Stacking_Masked_MasterGrid_.pdf'
name_dict['Stackname_nomask'] = 'Stacking_MasterGrid_.pdf'
name_dict['Average_Bin_Value'] = '_Average_R23_O32_Values.tbl'
name_dict['bin_fit_fits']='bin_emission_line_fit.fits'
name_dict['temp_fits_file']='_temperatures_metalicity.fits'
name_dict['temp_metallicity_pdf']='_Temp_Composite_Metallicity.pdf'



from Zcalbase_gal.Analysis.DEEP2_R23_O32 import Binning_and_Graphing_MasterGrid, Stackboth_MasterGrid, zoom_and_gauss_general, hstack_tables,  adaptivebinning, Stacking_voronoi, R_temp_calcul, calibration_plots, verification_tables
from Zcalbase_gal.Analysis.DEEP2_R23_O32.Plotting import more_plots, line_ratio_plotting, te_metal_plots

from Zcalbase_gal.Analysis import local_analog_calibration, green_peas_calibration

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf, intersect

#fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/n_split/'

from Metallicity_Stack_Commons.Metallicity_Stack_Commons import exclude_outliers, dir_date,lambda0, line_type, line_name, valid_table
from Metallicity_Stack_Commons.Metallicity_Stack_Commons.column_names import filename_dict
from Metallicity_Stack_Commons.Metallicity_Stack_Commons.plotting import balmer
from Metallicity_Stack_Commons.Metallicity_Stack_Commons.analysis import attenuation, composite_indv_detect, error_prop






'''outfile01 = fitspath+ 'Arrays_R23O32bin01MasterGrid.npz'
    outfile025 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz' #this file has the average R23 and O32 values for grid method
    outsingle_O32 = fitspath +'single_grid_O32.npz'
    outsingle_R23 = fitspath +'single_grid_R23.npz'
    outdouble_bin = fitspath +'nsplit_grid.npz'
    if dataset =='Grid' : grid_data_file = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'     ##np.load(outfile025)   ###This will have to be changed if we start doing the 01 analysis again (but we haven't worked on that analysis in a year) 
    if dataset == 'O32_Grid': grid_data_file = fitspath +'single_grid_O32.npz'      # = np.load(outsingle_O32)
    if dataset == 'R23_Grid': grid_data_file = fitspath +'single_grid_R23.npz'      # = np.load(outsingle_R23)
    #else: grid_data = np.load(outfile01)
    #if dataset == 'Double_Bin': grid_data_file = fitspath +'double_grid.npz'''




#############Getting username##############

import getpass

username = getpass.getuser()
print(username)
if username == 'reagenleimbach':
    fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'
if username == 'fill in':
    fitspath_ini = 'fill in '


def get_det3(fitspath, individual_detect=False):
    for ii in range(1,5):
        file1 = fitspath_ini+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

    objno = data0['OBJNO']

    #Excluding Outliers 
    exclude_flag = exclude_outliers(objno)
    print("exclude flag: ", np.where(exclude_flag == 1))

    O2_ini = data0['OII_FLUX_MOD'].data
    O3_ini = 1.33*data0['OIIIR_FLUX_MOD'].data
    O4959_ini = data0['OIIIB_FLUX_MOD'].data
    O5007_ini = data0['OIIIR_FLUX_MOD'].data
    Hgamma_ini = data0['HG_FLUX_MOD'].data
    O4363_ini = data0['OIIIA_FLUX_MOD'].data
    Hdelta_ini = data0['HD_FLUX_MOD'].data
    
    Hb_ini = data0['HB_FLUX_MOD'].data
    R23_ini = (O2_ini+O3_ini)/Hb_ini
    O32_ini = O3_ini/O2_ini

    lR23_ini =  np.log10(R23_ini)
    lO32_ini = np.log10(O32_ini)

    SNR2_ini = data0['OII_SNR'].data
    SNR3_ini = data0['OIIIR_SNR'].data
    SNRH_ini = data0['HB_SNR'].data
    SNRHg_ini = data0['HG_SNR'].data
    SNR4363_ini = data0['OIIIA_SNR'].data

    
    print('O2 len:', len(O2_ini))

    #################################################################################
    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini>0) & (exclude_flag==0) & (lR23_ini < 1.4))[0]  #May limit the logR23 value further to 1.2, check the velocity dispersions of the high R23 spectra

    # Organize the R23_032 data
    data3 = data0[det3]

    R23 = R23_ini[det3]
    O32 = O32_ini[det3]
    lR23 = lR23_ini[det3]
    lO32 = lO32_ini[det3]
    
    Hb  = Hb_ini[det3]
    O2  = O2_ini[det3]
    O3  = O3_ini[det3]
    Hgamma = Hgamma_ini[det3]
    Hdelta = Hdelta_ini[det3]
    O4363 = O4363_ini[det3]
    O4959 = O4959_ini[det3]
    O5007 = O5007_ini[det3]
    SNR2 =SNR2_ini[det3]
    SNR3 =SNR3_ini[det3]
    SNRH =SNRH_ini[det3]
    SNRHG = SNRHg_ini[det3]
    SNR4363 = SNR4363_ini[det3]
    individual_names = objno[det3]

    

    n2= ('ID','logR23','logO32','OII_3727_Flux_Gaussian','O3_Flux_Gaussian','HGAMMA_Flux_Gaussian','HDELTA_Flux_Gaussian',
         'OIII_4363_Flux_Gaussian','OIII_4958_Flux_Gaussian','OIII_5007_Flux_Gaussian',
         'HBETA_Flux_Gaussian', 'O2_S/N', 'O3_S/N', 'RH_S/N', 'HGAMMA_S/N', 'O4363_S/N')
    tab1 = Table([individual_names, lR23, lO32, O2, O3, Hgamma, Hdelta, O4363, O4959, O5007,
                  Hb, SNR2, SNR3, SNRH, SNRHG,SNR4363], names=n2)

    # We can create two different kinds of tables here of the R23_032 data (det3)
    asc.write(tab1, fitspath+'individual_properties.tbl', format='fixed_width_two_line')  #used to be get_det3_table.tbl
    #tab1.write(fitspath_ini+'get_det3_table.fit', format = 'fits', overwrite = True)
    

    # Saving the detection condition (np.where) statements that will be used to call the data
    if individual_detect == True:
        DEEP2_m = asc.read(fitspath_ini+'results_deeps_revised.tbl')
        temp_x = DEEP2_m['best.stellar.m_star'].data
        ML_data_det = np.where((np.isfinite(temp_x)==True)
                       & (temp_x>0) & (temp_x < 1e13) & (exclude_flag == 0))[0]

        #define column names
        n_comb = ('Mass_Luminosity_Individual', 'R23_O32_Individual')

        ###Making the columns the same lengths by adding np.nan to shorter dataset
        if len(ML_data_det) > len(det3):
            tot_nan = len(ML_data_det)-len(det3)
            nan_array = np.empty(tot_nan)
            nan_array[:] = np.nan
            new_det3 = np.concatenate((det3,nan_array))
            det_tab = Table([ML_data_det,new_det3], names=n_comb)
        else:
            tot_nan = len(det3)-len(ML_data_det)
            nan_array = np.empty(tot_nan)
            nan_array[:] = np.nan
            new_det = np.concatenate((ML_data_det,nan_array))
            det_tab = Table([new_det,det3], names=n_comb)
            
        




        print('Using mass evolve table from 09/24')
        mass_valid_idx_tab = np.load('/Users/reagenleimbach/Desktop/Zcalbase_gal/From_Caroline/mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz')
        mass_LHbeta_file = np.load('/Users/reagenleimbach/Desktop/Zcalbase_gal/From_Caroline/revised_75_112_113_300_600_1444_1444_mass_SFR_data.npz')
        
        bin_idx = mass_valid_idx_tab['bin_ind']
        mass= mass_LHbeta_file['mass']
        LHbeta = mass_LHbeta_file['lum']
        
        valid_idx = []
        for i in range(len(bin_idx)):           #bin_ind is 2D array of indices (an array for each bin)
            valid_idx += list(bin_idx[i])               #creates 1D array of all indices
        idx_array = np.array(valid_idx)
            
        mass_array =np.log10(mass[valid_idx])
        LHbeta_array = LHbeta[valid_idx]
        objectname = objno[valid_idx]
        

        n_massL = ('Object', 'Mass', 'Luminosity')
        mL_tab = Table([objectname, mass_array, LHbeta_array], names=n_massL)

        asc.write(det_tab,fitspath+'indexed_individual.tbl',format = 'fixed_width_two_line')
        asc.write(mL_tab, fitspath+'mass_LHbeta_tab.tbl', format = 'fixed_width_two_line')
     




    return individual_names, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3    #, O2_det3, O3_det3, Hb_det3



#redue calling for look at hb

def gettime(org_name,fitspath_ini):
    today = date.today()
    path0 = org_name+'_' + "%02i%02i" % (today.month, today.day)+'/'
    if not exists(path0):
        os.mkdir(fitspath_ini+path0)
        fitspath= fitspath_ini+path0
    else: print("Path already exists")
    print(fitspath)
    return fitspath

def check_verification_table(fitspath_ini, dataset, combine_flux_ascii):
    verification_table = fitspath_ini+'verification_tables/'+dataset+'_verification_tab.tbl'
    if exists(verification_table):
        return verification_table
    else:
        print('Making verification table')
        verification_tables.verification_master(fitspath,dataset, combine_flux_ascii)
        return verification_table


def run_grid_R23_O32_analysis(dataset,y_correction, n_split, adaptive = False, dustatten = 'False', mask='None'):
    #dataset options: Grid, O32_Grid, R23_Grid, n_Bins

    if dataset == 'O32_Grid': org_name = 'O32_Grid'
    if dataset == 'R23_Grid': org_name = 'R23_Grid'
    if dataset == 'Grid': org_name = 'R23O32_Grid'
    if dataset == 'n_Bins': org_name = 'R23O32_Manual'
    fitspath = gettime(org_name,fitspath_ini)
    
    individual_ID, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3(fitspath)     #, O2_det3, O3_det3, Hb_det3

    print("length R23: ", len(R23))
    

    #Each bin will be split in half
    #Must sum to 2799 for Zcalbase_gal Analysis
    if adaptive == False: galinbin = [400,400,400,400,400,400,409] 
    if adaptive == True: galinbin = [458,450,400,300,300,275,250,200,176] 
    print('# of Gal in Bin:', galinbin)

    Bin_pdf_pages = fitspath +dataset+ name_dict['gridpdf']
    Bin_outfile = fitspath +dataset+name_dict['gridnpz']
    
    
    if dataset =='O32_Grid': 
        Binning_and_Graphing_MasterGrid.single_grid_O32(fitspath, Bin_pdf_pages, Bin_outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, galinbin, adaptive)    

    if dataset =='R23_Grid':
        Binning_and_Graphing_MasterGrid.single_grid_R23(fitspath, Bin_pdf_pages, Bin_outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin)    


    if dataset =='Grid': 
        R23_bin = 0.25
        O32_bin = 0.25
        Binning_and_Graphing_MasterGrid.making_Grid(fitspath, Bin_pdf_pages,Bin_outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, R23_bin, O32_bin)   



    if dataset == 'n_Bins':
        pdf_pages = fitspath +'n_Bins_grid.pdf'
        grid_data_file = fitspath +'n_Bins_grid.npz'
        Binning_and_Graphing_MasterGrid.n_times_binned(fitspath,
                                                       Bin_pdf_pages,
                                                       Bin_outfile,
                                                       n_split,
                                                       individual_ID,
                                                       R23,
                                                       O32,
                                                       O2,
                                                       O3,
                                                       Hb,
                                                       SNR2,
                                                       SNR3,
                                                       SNRH,
                                                       det3,
                                                       data3,
                                                       galinbin,
                                                       adaptive)


    print('made npz, pdf files , testmastergrid(need to find if this is used anywhere)')
    print('finished Binning_and_Graphing_MasterGrid')

    



    #Stackboth_MasterGrid
    #Option to Change: Bin size  
    

    #Option to Change: Masking the night sky emission lines 
    if dataset == 'Grid': 
        if mask == True:
            Stack_name = dataset+ name_dict['Stackname']
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,Bin_outfile)
        else:
            Stack_name = dataset+name_dict['Stackname_nomask']
            Stackboth_MasterGrid.run_Stacking_Master(fitspath,Stack_name,Bin_outfile)
    else:
        if mask == True:
            Stack_name = dataset+ name_dict['Stackname']
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,Bin_outfile)
        else:
            Stack_name = dataset+name_dict['Stackname_nomask']
            Stackboth_MasterGrid.run_Stacking_Master(fitspath, Stack_name,Bin_outfile)

    #Outfile and pdf both use name
    print('finished with stacking,' + Stack_name + 'pdf and fits files created')



    #Zoom_and_gauss_general
    
    outfile_grid = fitspath + filename_dict['comp_spec']
    print(outfile_grid)
    stack2D, header = fits.getdata(outfile_grid, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    #Spect_1D = fits.getdata(outfile_grid)
    dispersion = header['CDELT1']
    binning_avg_asc = fitspath+ filename_dict['bin_info']  

    indv_bin_info = fitspath+ filename_dict['indv_bin_info']  #used to be 'n_Bins_2d_binning_datadet3.tbl'
    
    lineflag = np.zeros(len(wave))
    for ii in lambda0:   
        idx = np.where(np.absolute(wave - ii)<=5)[0]
        if len(idx) > 0:
            lineflag[idx] = 1
    #Option to change: Constants used as initial guesses for gaussian fit
    
    s= 1.0
    a= 1.0
    c= 2.0
    s1= 1.3
    a1= 1.5
    s2= 5.0
    a2= 1.8
    

    zoom_and_gauss_general.zm_general(dataset, fitspath,stack2D, wave, lineflag, dispersion, y_correction, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)

    print('finished gaussian fitting:,' +fitspath+'_'+dataset+'_Zoomed_Gauss_* pdfs and fits created')

    #hstack_table
    #Option to change: name of new fits file created

    intro = fitspath + '/' + dataset + name_dict['Average_Bin_Value']
    asc_intro = asc.read(intro)
    table_files = glob.glob(fitspath +'/' + dataset + '_flux_gaussian_*.tbl') 
    combine_flux_table = fitspath + name_dict['bin_fit_fits']
    combine_flux_ascii = fitspath + filename_dict['bin_fit']

    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)
        
    print('combine_flux_table created')

        
            
    ########FIX THIS CODE##########line_ratio_plotting
    #I need to go back through and figure out what is the average and what is the composite
    #line_ratio_plotting.Plotting_Data1(fitspath,dataset,combine_flux_ascii, binning_avg_asc)

    #Verification Table 
    valid_table.make_validation_table(fitspath)
    verification_table = fitspath + filename_dict['bin_valid']
    
    valid_table.compare_to_by_eye(fitspath,dataset)
    verification_table_revised = fitspath + filename_dict['bin_valid_rev']
    
    #R_temp_calcul
    #Not going to run the R_temp_calcul.run_function for the 'bin_valid' table because these values
    #are proven to be incomplete. The bin_valid_rev table is different. 
    R_temp_calcul.run_function(fitspath, dataset,verification_table_revised, dustatt= False)

    if dustatten == True :
        balmer.HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits', use_revised=False) # outfile_grid,

        attenuation.EBV_table_update(fitspath, use_revised= False)

        #non_atten_value_table = asc.read(temp_m_gascii)
        #EBV_HgHb = non_atten_value_table['EBV_HgHb']
        #temp_tab_revised = fitspath+ filename_dict['bin_derived_prop_rev']
        R_temp_calcul.run_function(fitspath, dataset, verification_table_revised, dustatt= True)
        
    '''
    #####Check Dust Attenuation#####
    temp = fitspath + filename_dict['bin_derived_prop']
    temp_revised = fitspath +filename_dict['bin_derived_prop_dustcorr']

    ttab = asc.read(temp)
    trtab = asc.read(temp_revised)
    
    metallicity = ttab['12+log(O/H)'].data
    metallicity_r = trtab['12+log(O/H)'].data

    print('Metallicity' , metallicity)
    print('##################################################################')
    print('Dust Attenuated Metallicity', metallicity_r)'''

    ####Error Propagation#####
    error_prop.fluxes_derived_prop(fitspath, binned_data = True, revised = True)

    
    ####Individual Detections#####
    #composite_indv_detect.main(fitspath, dataset= '', revised = False, det3=True)
    #print('Individual Detections Complete')

    ###Te_metal Plots####
    
    #te_metal_plots.plotting_te_metal(fitspath, revised=False)
    #te_metal_plots.plotting_te_metal(fitspath,  revised=True)


    ###Calibration Plots###
    #calibration_plots.LAC_GPC_plots(fitspath, dataset, revised= False)
    calibration_plots.LAC_GPC_plots(fitspath, dataset, revised= True)
    
'''

    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    #asc_table_det3 = asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    
    temp_m_gascii = fitspath+ 'n_Bins_temperatures_metalicity.tbl'
    
    
    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_m_gascii, verification_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_m_gascii, verification_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_gascii, verification_table)
    more_plots.hist_for_bin(dataset, asc_table2)
    '''






def run_voronoi_R23_O32_analysis(dataset,y_correction, mask='None'):
    #dataset options: Voronoi10,Voronoi14, Voronoi20
    
    if dataset == 'Voronoi10': org_name = 'Voronoi10'
    if dataset == 'Voronoi14': org_name = 'Voronoi14'
    if dataset == 'Voronoi20': org_name = 'Voronoi20'
    if dataset == 'Double_Bin': org_name = 'Double_Bin'
    fitspath = gettime(org_name,fitspath_ini)
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3(fitspath) # , O2_det3, O3_det3, Hb_det3
    
    #Grid Methods of Organizing Data
    #Both One and Two Dimensional Analyzes 
    #Stackboth_MasterGrid
    #Option to Change: Bin size 


    ###Voronoi###


    #Adaptive Binning 
    #Options to Change: Signal/Noise Size
    if dataset == 'Voronoi10':
        sn_size = 10
        txt_file = fitspath + 'voronoi10_2d_binning_outputs.txt'
        asc_table1 = fitspath+'bin_info.tbl' #used to be called asc_tab_fill in name
        asc_table2 = fitspath+'voronoi10_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_10.tbl'
        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)   #, O2_det3, O3_det3, Hb_det3)
        
    if dataset == 'Voronoi14':
        sn_size = 14
        txt_file = fitspath + 'voronoi14_2d_binning_outputs.txt'
        asc_table1 = fitspath+'bin_info.tbl'
        asc_table2 = fitspath+'voronoi14_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_14.tbl'
        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)    #, O2_det3, O3_det3, Hb_det3
        
    if dataset == 'Voronoi20':
        sn_size = 20
        txt_file = fitspath + 'voronoi20_2d_binning_outputs.txt'
        asc_table1 = fitspath+'bin_info.tbl'
        asc_table2 = fitspath+'voronoi20_2d_binning_datadet3.tbl' 

        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)   #, O2_det3, O3_det3, Hb_det3

    if dataset == 'Double_Bin':
        galinbin = 400
        pdf_pages = fitspath +'double_grid.pdf'
        outfile = fitspath +'double_grid.npz'
        asc_table1 = fitspath+ '/bin_info.tbl'
        asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
        Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin)   #, O2_det3, O3_det3, Hb_det3

    
    #Stacking_voronoi
    #Option to Change: 
    if dataset == "Voronoi10": voronoi_data = fitspath+ 'voronoi10_2d_binning_datadet3.tbl'
    if dataset == "Voronoi14": voronoi_data = fitspath + 'voronoi14_2d_binning_datadet3.tbl'
    if dataset == "Voronoi20": voronoi_data = fitspath + 'voronoi20_2d_binning_datadet3.tbl'

    if dataset == 'Double_Bin': grid_data = np.load(outfile)
    if dataset == 'Double_Bin': voronoi_data = asc_table2

    print('### outfile for datadet3: '+voronoi_data)
        
    #Option to Change: Masking the night sky emission lines
    ######Check to make sure tables are going into right places#####
    if mask == True:
        Stack_name = 'Stacking'+dataset+'_output.pdf'
        Stacking_voronoi.run_Stacking_Master_mask(fitspath_ini, dataset, fitspath, voronoi_data, det3, asc_table1, Stack_name)
    else:
        Stack_name = 'Stacking'+dataset+'_output.pdf'
        Stacking_voronoi.run_Stacking_Master(fitspath_ini, dataset, fitspath, voronoi_data, det3, Stack_name)

    #Outfile and pdf both use name
    print('finished with stacking,' + Stack_name + ' pdf and fits files created')



    #Zoom_and_gauss_general

    Stack_name= Stack_name.replace('.pdf', '.fits')
    outfile_voronoi = fitspath+ Stack_name
    stack2D, header = fits.getdata(outfile_voronoi,header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    #Spect_1D = fits.getdata(outfile_voronoi)
    dispersion = header['CDELT1']

    lineflag = np.zeros(len(wave))
    for ii in lambda0:   
        idx = np.where(np.absolute(wave - ii)<=5)[0]
        if len(idx) > 0:
            lineflag[idx] = 1

    ###Delete after you make sure it is not used###
    '''#tab = asc_table1
    if dataset == 'Voronoi10':
        tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_14.tbl'
       #asc_tab = asc.read(tab)
    else: 
        tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_10.tbl'
        #asc_tab = asc.read(tab)'''
            
    #Option to change: Constants used as initial guesses for gaussian fit
    s= 1.0
    a= 1.0
    c= 2.0
    s1= 1.3
    a1= 4.7
    s2= 10.0
    a2= -2.0
    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave,lineflag, dispersion, y_correction, s,a,c,s1,a1,s2,a2,tab=asc_table1)

    print('finished gaussian emission fitting pdf and tables created')

    #hstack_table
    #Option to change: name of new fits file created
    intro = fitspath + dataset+'_Average_R23_O32_Values.tbl'
    asc_intro = asc.read(intro)
    table_files = glob.glob(fitspath+ dataset+'_flux_gaussian_*.tbl')
    combine_flux_fits = fitspath+'bin_emission_line_fit.fits'
    combine_flux_ascii = fitspath+'bin_emission_line_fit.tbl'
    print(combine_flux_ascii)
    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)

    print(dataset+'_combine_flux_table created')

    ####### FIX THIS PLOTS ######line_ratio_plotting
    #I need to go back through and figure out what is the average and what is the composite
    line_ratio_plotting.Plotting_Data1(fitspath, dataset, combine_flux_ascii, asc_table1)


        
    #R_temp_calcul
        
    temp_met_ascii = fitspath+ '/'+dataset+'_temperatures_metalicity.tbl'
    temp_met_fits = fitspath+ '/'+dataset+'_temperatures_metalicity.fits'
    pdf_name_temp_met = dataset+'_Temp_Composite_Metallicity.pdf'

    R_temp_calcul.run_function(fitspath, dataset, temp_m_gascii , temp_m_gfits, temp_m_gpdf_name, combine_flux_ascii, dustatt=False)

  
    ###Calibration Plots###
    calibration_plots.LAC_GPC_plots(fitspath,dataset,temp_met_ascii)

    '''
    ###Verification Tables###
    ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    verification_tables.verification(fitspath, dataset, temp_met_ascii, combine_flux_ascii, ver_tab)'''

    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    #asc_table_det3 = asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    m_ver_table = fitspath_ini+'Master_Verification_Tables/'+dataset+'_table.tbl'
    #ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    
    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.hist_for_bin(dataset, asc_table2)
    
    












def run_two_times_binned_analysis(dataset,y_correction, adaptive = False, mask='None'):
    #dataset must equal Double_bin
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3()


    if dataset == 'Double_Bin':
        if adaptive == False: galinbin = [400,400,400,400,400,400,409] #Each bin will be split in half
        if adaptive == True: galinbin = [458,450,400,300,300,275,250,200,176] #Must sum to 2809 
        pdf_pages = fitspath +'double_grid.pdf'
        grid_data_file = fitspath +'double_grid.npz'
        asc_table1 = fitspath+ '/bin_info.tbl'
        asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
        Binning_and_Graphing_MasterGrid.two_times_binned(fitspath,
                                                         pdf_pages,
                                                         grid_data_file,
                                                         R23,
                                                         O32,
                                                         O2,
                                                         O3,
                                                         Hb,
                                                         SNR2,
                                                         SNR3,
                                                         SNRH,
                                                         det3,
                                                         data3,
                                                         galinbin,
                                                         adaptive) 
        

    #Stacking_MASKED_MASTERGRID
    Stack_name = 'Stacking_Masked_MasterGrid_single'+dataset+'.pdf' 
    Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data_file)

    #Outfile and pdf both use name
    print('finished with stacking,' + Stack_name + ' pdf and fits files created')



    #Zoom_and_gauss_general

    Stack_name= Stack_name.replace('.pdf', '.fits')
    outfile_voronoi = fitspath+ Stack_name
    stack2D, header = fits.getdata(outfile_voronoi,header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    #Spect_1D = fits.getdata(outfile_voronoi)
    dispersion = header['CDELT1']

    lineflag = np.zeros(len(wave))
    for ii in lambda0:   
        idx = np.where(np.absolute(wave - ii)<=5)[0]
        if len(idx) > 0:
            lineflag[idx] = 1

    ###Delete after you make sure it is not used###
    #tab = asc_table1
    if dataset == 'Voronoi10':
        tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_14.tbl'
       #asc_tab = asc.read(tab)
    else: 
        tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_10.tbl'
        #asc_tab = asc.read(tab)
            
    #Option to change: Constants used as initial guesses for gaussian fit
    s=1.0
    a= 1.0
    c = 2.0
        
    s1=1.3
    a1= 4.7
    s2 = 10.0
    a2 = -2.0
    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave,lineflag, dispersion, y_correction, s,a,c,s1,a1,s2,a2,tab=asc_table1)

    print('finished gaussian emission fitting pdf and tables created')

    #hstack_table
    #Option to change: name of new fits file created
    intro = fitspath + dataset+'_Average_R23_O32_Values.tbl'
    asc_intro = asc.read(intro)
    table_files = glob.glob(fitspath+ dataset+'_flux_gaussian_*.tbl')
    combine_flux_fits = fitspath+'bin_emission_line_fit.fits'
    combine_flux_ascii = fitspath+'bin_emission_line_fit.tbl'
    print(combine_flux_ascii)
    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)

    print(dataset+'_combine_flux_table created')

    ####### FIX THIS PLOTS ######line_ratio_plotting
    #I need to go back through and figure out what is the average and what is the composite
    line_ratio_plotting.Plotting_Data1(fitspath, dataset, combine_flux_ascii, asc_table1)


        
    #R_temp_calcul
        
    temp_met_ascii = fitspath+ '/'+dataset+'_temperatures_metalicity.tbl'
    temp_met_fits = fitspath+ '/'+dataset+'_temperatures_metalicity.fits'
    pdf_name_temp_met = dataset+'_Temp_Composite_Metallicity.pdf'

    R_temp_calcul.run_function(fitspath,dataset, temp_met_ascii,temp_met_fits, pdf_name_temp_met, combine_flux_ascii)
  
    ###Calibration Plots###
    calibration_plots.LAC_GPC_plots(fitspath,dataset,temp_met_ascii)

    '''
    ###Verification Tables###
    ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    verification_tables.verification(fitspath, dataset, temp_met_ascii, combine_flux_ascii, ver_tab)'''

    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    #asc_table_det3 = asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    m_ver_table = fitspath_ini+'Master_Verification_Tables/'+dataset+'_table.tbl'
    #ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    
    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_ascii, m_ver_table)
    more_plots.hist_for_bin(dataset, asc_table2)











#Below function will run the individual functions in the codes above that produce graphs
#Enter a keyword for want to indictate what function you want to run
#This will ease the reproduction process
#CHECK: function defaults to put new graphs in fitspath. Make sure you don't over write something you need
def run_individual_functions(fitspath, want, adaptive, y_correction, dustatten= False, individual=False):
    #Keywords: binning_and_graphing, stack_mastergrid, zoom, R_cal_temp, line_ratio_plotting
    
    dataset = 'n_Bins'
    


    if want == 'binning_and_graphing':
        R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3(fitspath)
        if adaptive == False: galinbin = [400,400,400,400,400,400,409] #Each bin will be split in half
        if adaptive == True: galinbin = [458,450,400,300,300,275,250,200,176] #Must sum to 2800 
        pdf_pages = fitspath +'n_Bins_grid.pdf'
        grid_data_file = fitspath +'n_Bins_grid.npz'
        asc_table1 = fitspath+ '/bin_info.tbl'
        asc_table2 = fitspath+ 'n_Bins_2d_binning_datadet3.tbl'
        Binning_and_Graphing_MasterGrid.n_times_binned(fitspath,
                                                       pdf_pages,
                                                       grid_data_file,
                                                       n_split,
                                                       R23,
                                                       O32,
                                                       O2,
                                                       O3,
                                                       Hb,
                                                       SNR2,
                                                       SNR3,
                                                       SNRH,
                                                       det3,
                                                       data3,
                                                       galinbin,
                                                       adaptive)



    if want == 'stack_mastergrid':
        outdouble_bin = fitspath +'nsplit_grid.npz'
        grid_data_file = fitspath +'n_Bins_grid.npz'
        Stack_name = 'Stacking_Masked_MasterGrid_'+dataset+'.pdf'
        Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data_file)



    if want == 'zoom':
        Stack_name = 'Stacking_Masked_MasterGrid_'+dataset+'.fits'
        outfile_grid = fitspath + Stack_name
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        dispersion = header['CDELT1']
        binning_avg_asc = fitspath+'/bin_info.tbl'
    
    
        lineflag = np.zeros(len(wave))
        for ii in lambda0:   
            idx = np.where(np.absolute(wave - ii)<=5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1
    
    
        s= 1.0
        a= 1.0
        c= 2.0
        s1= 1.3
        a1= 1.5
        s2= 5.0
        a2= 1.8
    

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, wave, lineflag, dispersion, y_correction, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)


    if want == 'R_cal_temp':
        combine_flux_ascii = fitspath + 'bin_emission_line_fit.tbl'
        temp_m_gascii = fitspath+ '/nsplit_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/nsplit_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'nsplit_Temp_Composite_Metallicity.pdf'

        if dustatten == 'False': R_temp_calcul.run_function(fitspath, dataset, temp_m_gascii , temp_m_gfits, temp_m_gpdf_name, combine_flux_ascii, dustatt= False)
        if dustatten == 'True': R_temp_calcul.run_function(fitspath, dataset, temp_m_gascii , temp_m_gfits, temp_m_gpdf_name, combine_flux_ascii, dustatt= True)

        
    if want =='line_ratio_plotting': 
        combine_flux_ascii = fitspath + 'bin_emission_line_fit.tbl'
        binning_avg_asc = fitspath+'/bin_info.tbl'
        line_ratio_plotting.Plotting_Data1(fitspath,dataset,combine_flux_ascii, binning_avg_asc)


    if want =='calibration_plots':
        temp_m_gascii = fitspath+ '/nsplit_temperatures_metalicity.tbl'
        calibration_plots.LAC_GPC_plots(fitspath, dataset, temp_m_gascii)


    print(want, 'is done')

    










#####Not Using Yet######
'''if individual == 'True':
            individual_variables_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0902/Individual_ratio_temperature.tbl'
            R_temp_calcul.run_individual_function(fitspath,dataset,individual_variables_ascii, 
'''
