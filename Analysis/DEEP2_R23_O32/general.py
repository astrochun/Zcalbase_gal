#This function runs the entire process start to finish
#This weekend combine the grid and voronoi if statements, voronoi20, log plots for met
#EW values:equival width
### 

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

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import Binning_and_Graphing_MasterGrid, Stackboth_MasterGrid, zoom_and_gauss_general, hstack_tables,  adaptivebinning, Stacking_voronoi, R_temp_calcul, line_ratio_plotting,calibration_plots, verification_tables, more_plots

from Zcalbase_gal.Analysis import local_analog_calibration, green_peas_calibration

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Double_Bin_Grid_Stacking_0722/'
fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]

'''Complete List
lambda0 =[3726.16, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  #, '4101.73, 4363.21, 4861.32']
        #[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single', 'Balmer', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','H_10', 'H_9','NeIII','Hzeta','HeI','3967', 'HEPSIL', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  '''

lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  #, '4101.73, 4363.21, 4861.32']
        #[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Oxy2', 'Single','Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  

def exclude_outliers(objno):
    flag = np.zeros(len(objno), dtype=int)
    bad_data = np.array(['32007727', '32101412', '42006031', '32035286', '14023705'])
    for ii in range(len(bad_data)):
        idx = [xx for xx in range(len(objno)) if bad_data[ii] == str(objno[xx])][0]
        flag[idx] = 1
    
    return flag

def get_det3():
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
    print "exclude flag: ", np.where(exclude_flag == 1)
    
    O2_ini = data0['OII_FLUX_MOD']
    O3_ini = 1.33*data0['OIIIR_FLUX_MOD']
    Hb_ini = data0['HB_FLUX_MOD']
    R23_ini = (O2_ini+O3_ini)/Hb_ini
    O32_ini = O3_ini/O2_ini

    SNR2_ini = data0['OII_SNR']
    SNR3_ini = data0['OIIIR_SNR']
    SNRH_ini = data0['HB_SNR']

    logR23 = np.log10(R23)
    print 'O2 len:', len(O2_ini)
    
    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini>0) & (exclude_flag==0) & (logR23 < 1.4))[0]

    
    data3 = data0[det3]

    R23 = R23_ini[det3]
    O32 = O32_ini[det3]
    Hb  = Hb_ini[det3]
    O2  = O2_ini[det3]
    O3  = O3_ini[det3]
    SNR2 =SNR2_ini[det3]
    SNR3 =SNR3_ini[det3]
    SNRH =SNRH_ini[det3]

    print 'R23:',len(R23)
    O2_det3 =O2_ini[det3]
    O3_det3 =O3_ini[det3]
    Hb_det3 =Hb_ini[det3]

    n2= ('R23','O32', 'O2', 'O3', 'Hb', 'SNR2', 'SNR3', 'SNRH', 'det3', 'O2_det3', 'O3_det3', 'Hb_det3')  #'data3',
    tab1 = Table([R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, O2_det3, O3_det3, Hb_det3], names=n2)  #data3
    asc.write(tab1, fitspath+'get_det3_table2.tbl', format='fixed_width_two_line')
    #tab1.write(fitspath_ini+'get_det3_table.fit', format = 'fits', overwrite = True)
    
    return R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3    #, O2_det3, O3_det3, Hb_det3



#redue calling for look at hb


def run_grid_R23_O32_analysis(dataset,y_correction, adaptive = False, mask='None'):
    #dataset options: Grid, O32_Grid, R23_Grid, Double_Bin
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3()     #, O2_det3, O3_det3, Hb_det3
    
    #Grid Methods of Organizing Data
    #Both One and Two Dimensional Analyzes 
    #Stackboth_MasterGrid
    #Option to Change: Bin size 

    
    #Binning and Graphing MasterGrid
    #Options to Change: Bin Size
    #galinbin = 400
    if adaptive == False: galinbin = [400,400,400,400,400,400,409] #Each will be split in half
    if adaptive == True: galinbin = [908,700,575,450,176] #Must sum to 2809 
    print '# of Gal in Bin:', galinbin
    if dataset =='O32_Grid': 
        pdf_pages = fitspath +'single_grid_O32.pdf'
        outfile = fitspath +'single_grid_O32.npz'
        Binning_and_Graphing_MasterGrid.single_grid_O32(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, galinbin, adaptive)    #, O2_det3, O3_det3, Hb_det3,

    if dataset =='R23_Grid':
        pdf_pages = fitspath +'single_grid_R23.pdf'
        outfile = fitspath +'single_grid_R23.npz'
        Binning_and_Graphing_MasterGrid.single_grid_R23(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin)     #O2_det3, O3_det3, Hb_det3,


    if dataset =='Grid': 
        pdf_pages_hex = PdfPages(fitspath+'R23_O32_bin025_scatter_and_hexbin_MasterGrid.pdf')
        outfile1_npz = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'
        R23_bin = 0.25
        O32_bin = 0.25
        binstr = 025
        
        Binning_and_Graphing_MasterGrid.making_Grid(fitspath, pdf_pages_hex,outfile1_npz,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, R23_bin, O32_bin)   #O2_det3, O3_det3, Hb_det3,



    if dataset == 'Double_Bin':
        if adaptive == False: galinbin = [400,400,400,400,400,400,409] #Each bin will be split in half
        if adaptive == True: galinbin = [458,450,400,300,300,275,250,200,176] #Must sum to 2809 
        pdf_pages = fitspath +'double_grid.pdf'
        grid_data_file = fitspath +'double_grid.npz'
        asc_table1 = fitspath+ '/Double_Bin_binning_averages.tbl'
        asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
        if n_split == 2: Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, grid_data_file,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive) 
        if n_split == 3: Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, n_split, grid_data_file,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive)


    print 'made npz, pdf files , testmastergrid(need to find if this is used anywhere)'
    print 'finished Binning_and_Graphing_MasterGrid'

    



    #Stackboth_MasterGrid
    #Option to Change: Bin size  
    outfile01 = fitspath+ 'Arrays_R23O32bin01MasterGrid.npz'
    outfile025 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz' #this file has the average R23 and O32 values for grid method
    outsingle_O32 = fitspath +'single_grid_O32.npz'
    outsingle_R23 = fitspath +'single_grid_R23.npz'
    outdouble_bin = fitspath +'double_grid.npz'
    if dataset =='Grid' : grid_data_file = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'     ##np.load(outfile025)   ###This will have to be changed if we start doing the 01 analysis again (but we haven't worked on that analysis in a year) 
    if dataset == 'O32_Grid': grid_data_file = fitspath +'single_grid_O32.npz'      # = np.load(outsingle_O32)
    if dataset == 'R23_Grid': grid_data_file = fitspath +'single_grid_R23.npz'      # = np.load(outsingle_R23)
    #else: grid_data = np.load(outfile01)
    #if dataset == 'Double_Bin': grid_data_file = fitspath +'double_grid.npz'

        
        
    #Option to Change: Masking the night sky emission lines 
    if dataset == 'Grid': 
        if mask == True:
            Stack_name = 'Stacking_Masked_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data_file)
        else:
            Stack_name = 'Stacking_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master(fitspath, Stack_name,grid_data_file)
    else:
        if mask == True:
            Stack_name = 'Stacking_Masked_MasterGrid_single'+dataset+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data_file)
        else:
            Stack_name = 'Stacking_MasterGrid__single'+dataset+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master(fitspath, Stack_name,grid_data_file)

    #Outfile and pdf both use name
    print 'finished with stacking,' + Stack_name + 'pdf and fits files created'



    #Zoom_and_gauss_general
    Stack_name = Stack_name.replace('.pdf', '.fits')
    outfile_grid = fitspath + Stack_name
    print outfile_grid
    stack2D, header = fits.getdata(outfile_grid, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    #Spect_1D = fits.getdata(outfile_grid)
    dispersion = header['CDELT1']
    binning_avg_asc = fitspath+'/'+dataset+'binning_averages.tbl'
    
    
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
    

    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name,  y_correction, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)

    print 'finished gaussian fitting:,' +fitspath+'_'+dataset+'_Zoomed_Gauss_* pdfs and fits created'

    #hstack_table
    #Option to change: name of new fits file created

    if dataset == 'Grid':
        intro = fitspath + 'Grid_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/Grid_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'Grid_combined_flux_table.fits'
        combine_flux_ascii = fitspath + 'Grid_combined_flux_table.tbl'

    if dataset == 'O32_Grid':
        intro = fitspath + 'O32_Grid_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/O32_Grid_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'O32_Grid_combined_flux_table.fits'
        combine_flux_ascii = fitspath + 'O32_Grid_combined_flux_table.tbl'

    if dataset == 'R23_Grid':
        intro = fitspath + 'R23_Grid_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/R23_Grid_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'R23_Grid_combined_flux_table.fits'
        combine_flux_ascii = fitspath + 'R23_Grid_combined_flux_table.tbl'

    if dataset == 'Double_Bin':
        intro = fitspath + 'Double_Bin_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/Double_Bin_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'Double_Bin_combined_flux_table.fits'
        combine_flux_ascii = fitspath + 'Double_Bin_combined_flux_table.tbl'

    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)
        
    print 'combine_flux_table created'

        
            
    ########FIX THIS CODE##########line_ratio_plotting
    #I need to go back through and figure out what is the average and what is the composite
    line_ratio_plotting.Plotting_Data1(fitspath,dataset,combine_flux_ascii, binning_avg_asc)

       

    #R_temp_calcul
    if dataset == 'Grid':
        combine_flux_ascii = fitspath + 'Grid_combined_flux_table.tbl'
        temp_m_gascii = fitspath+ '/Grid_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/Grid_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'Grid_Temp_Composite_Metallicity.pdf'
    if dataset == 'O32_Grid':
        combine_flux_ascii = fitspath + 'O32_Grid_combined_flux_table.tbl'
        temp_m_gascii = fitspath+ '/O32_Grid_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/O32_Grid_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'O32_Grid_Temp_Composite_Metallicity.pdf'

    if dataset == 'R23_Grid':
        combine_flux_ascii = fitspath + 'R23_Grid_combined_flux_table.tbl'
        temp_m_gascii = fitspath+ '/R23_Grid_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/R23_Grid_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'R23_Grid_Temp_Composite_Metallicity.pdf'

    if dataset == 'Double_Bin':
        combine_flux_ascii = fitspath + 'Double_Bin_combined_flux_table.tbl'
        temp_m_gascii = fitspath+ '/Double_Bin_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/Double_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'Double_Bin_Temp_Composite_Metallicity.pdf'

    R_temp_calcul.run_function(fitspath, dataset, temp_m_gascii , temp_m_gfits, temp_m_gpdf_name, combine_flux_ascii) 

    ###Calibration Plots###
    calibration_plots.LAC_GPC_plots(fitspath, dataset, temp_m_gascii)
    
    ###Verification Tables###
    #ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    #verification_tables.verification(fitspath, dataset, temp_m_gascii, combine_flux_ascii, ver_tab)


    ###Making More Plots###
    #asc_table = combine_flux_ascii
    #temp_table = temp_met_ascii
    #asc_table_det3 = asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
    m_ver_table = fitspath_ini+'Master_Verification_Tables/'+dataset+'_table.tbl'
    #ver_tab = fitspath+'/'+dataset+'_verification.tbl'
    
    more_plots.ew_plot_R23(fitspath, combine_flux_ascii, temp_met_gascii, m_ver_table)
    more_plots.ew_plot_O32(fitspath, combine_flux_ascii, temp_met_gascii, m_ver_table)
    more_plots.R23_vs_O32_color(fitspath, combine_flux_ascii, temp_met_gascii, m_ver_table)
    more_plots.hist_for_bin(dataset, asc_table2)







def run_voronoi_R23_O32_analysis(dataset,y_correction, mask='None'):
    #dataset options: Voronoi10,Voronoi14, Voronoi20
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3() # , O2_det3, O3_det3, Hb_det3
    
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
        asc_table1 = fitspath+'voronoi10_binning_averages.tbl' #used to be called asc_tab_fill in name
        asc_table2 = fitspath+'voronoi10_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_10.tbl'
        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)   #, O2_det3, O3_det3, Hb_det3)
        
    if dataset == 'Voronoi14':
        sn_size = 14
        txt_file = fitspath + 'voronoi14_2d_binning_outputs.txt'
        asc_table1 = fitspath+'voronoi14_binning_averages.tbl'
        asc_table2 = fitspath+'voronoi14_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_14.tbl'
        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)    #, O2_det3, O3_det3, Hb_det3
        
    if dataset == 'Voronoi20':
        sn_size = 20
        txt_file = fitspath + 'voronoi20_2d_binning_outputs.txt'
        asc_table1 = fitspath+'voronoi20_binning_averages.tbl'
        asc_table2 = fitspath+'voronoi20_2d_binning_datadet3.tbl' 

        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3)   #, O2_det3, O3_det3, Hb_det3

    if dataset == 'Double_Bin':
        galinbin = 400
        pdf_pages = fitspath +'double_grid.pdf'
        outfile = fitspath +'double_grid.npz'
        asc_table1 = fitspath+ '/Double_Bin_binning_averages.tbl'
        asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
        Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, outfile,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin)   #, O2_det3, O3_det3, Hb_det3

    
    #Stacking_voronoi
    #Option to Change: 
    if dataset == "Voronoi10": voronoi_data = fitspath+ 'voronoi10_2d_binning_datadet3.tbl'
    if dataset == "Voronoi14": voronoi_data = fitspath + 'voronoi14_2d_binning_datadet3.tbl'
    if dataset == "Voronoi20": voronoi_data = fitspath + 'voronoi20_2d_binning_datadet3.tbl'

    if dataset == 'Double_Bin': grid_data = np.load(outfile)
    if dataset == 'Double_Bin': voronoi_data = asc_table2

    print '### outfile for datadet3: '+voronoi_data
        
    #Option to Change: Masking the night sky emission lines
    ######Check to make sure tables are going into right places#####
    if mask == True:
        Stack_name = 'Stacking'+dataset+'_output.pdf'
        Stacking_voronoi.run_Stacking_Master_mask(fitspath_ini, dataset, fitspath, voronoi_data, det3, asc_table1, Stack_name)
    else:
        Stack_name = 'Stacking'+dataset+'_output.pdf'
        Stacking_voronoi.run_Stacking_Master(fitspath_ini, dataset, fitspath, voronoi_data, det3, Stack_name)

    #Outfile and pdf both use name
    print 'finished with stacking,' + Stack_name + ' pdf and fits files created'



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
    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave,lineflag, dispersion, lambda0, line_type, line_name,  y_correction, s,a,c,s1,a1,s2,a2,tab=asc_table1)

    print 'finished gaussian emission fitting pdf and tables created'

    #hstack_table
    #Option to change: name of new fits file created
    intro = fitspath + dataset+'_Average_R23_O32_Values.tbl'
    asc_intro = asc.read(intro)
    table_files = glob.glob(fitspath+ dataset+'_flux_gaussian_*.tbl')
    combine_flux_fits = fitspath+dataset+'_combined_flux_table.fits'
    combine_flux_ascii = fitspath+dataset+'_combined_flux_table.tbl'
    print combine_flux_ascii
    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)

    print dataset+'_combine_flux_table created'

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
    
    










###USING TO TEST IF I AM WRITING THE DOUBLE BIN METHOD CORRECTLY###


def TEST_FUNCTION_TO_RUN_TWO_BINS(dataset,y_correction, n_split, adaptive = False, mask='None'):
    #dataset options: Voronoi10,Voronoi14, Voronoi20
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = get_det3()


    if dataset == 'Double_Bin':
        if adaptive == False: galinbin = [400,400,400,400,400,400,409] #Each bin will be split in half
        if adaptive == True: galinbin = [458,450,400,300,300,275,250,200,176] #Must sum to 2809 
        pdf_pages = fitspath +'double_grid.pdf'
        grid_data_file = fitspath +'double_grid.npz'
        asc_table1 = fitspath+ '/Double_Bin_binning_averages.tbl'
        asc_table2 = fitspath+ 'Double_Bin_2d_binning_datadet3.tbl'
        if n_split == 2: Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, grid_data_file,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive) 
        if n_split == 3: Binning_and_Graphing_MasterGrid.two_times_binned(fitspath, pdf_pages, n_split, grid_data_file,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive) 

    #Stacking_MASKED_MASTERGRID
    Stack_name = 'Stacking_Masked_MasterGrid_single'+dataset+'.pdf' 
    Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data_file)

    #Outfile and pdf both use name
    print 'finished with stacking,' + Stack_name + ' pdf and fits files created'



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
    zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave,lineflag, dispersion, lambda0, line_type, line_name, y_correction, s,a,c,s1,a1,s2,a2,tab=asc_table1)

    print 'finished gaussian emission fitting pdf and tables created'

    #hstack_table
    #Option to change: name of new fits file created
    intro = fitspath + dataset+'_Average_R23_O32_Values.tbl'
    asc_intro = asc.read(intro)
    table_files = glob.glob(fitspath+ dataset+'_flux_gaussian_*.tbl')
    combine_flux_fits = fitspath+dataset+'_combined_flux_table.fits'
    combine_flux_ascii = fitspath+dataset+'_combined_flux_table.tbl'
    print combine_flux_ascii
    hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)

    print dataset+'_combine_flux_table created'

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

























###NOT IN USE###
'''
    temp_table= asc.read(temp_met_ascii)
    SN_4363 = temp_table['S/N_4363']
    det_4363 = np.where((SN_4363>=3))[0]
    nan_det_4363 = np.where((SN_4363<3))[0]
    
    print 'Begin Local analog Calibration'

    ###Implimenting Local analog calibration###
    derived = asc.read(fitspath_ini +'DEEP2_R23_O32_derived.tbl')
    derived_MACT = asc.read(fitspath_ini +'MACT_R23_O32_derived.tbl')

    out_pdf = fitspath+ '/'+dataset+'_LAC.pdf'
    O32_all = temp_table['O32_Composite']
    R23_all = temp_table['R23_Composite']
    com_O_log = temp_table['com_O_log']  #This is the 12+log(OH) value
    ID = temp_table['ID']
    print O32_all

    det_O32 = [O32_all[det_4363]]
    det_R23 = [R23_all[det_4363]]
    det_OH  = [com_O_log[det_4363]]
    det_ID  = [ID[det_4363]]

    nandet_O32 = [O32_all[nan_det_4363]]
    nandet_R23 = [R23_all[nan_det_4363]]
    nandet_OH  = [com_O_log[nan_det_4363]]
    nandet_ID  = [ID[nan_det_4363]]
    
    lR23 = [det_R23,nandet_R23]
    lO32 = [det_O32,nandet_O32]
    OH   = [det_OH,nandet_OH]

    local_analog_calibration.main(lR23, lO23, OH, out_pdf,ctype=['b','g'],label=['Detection','Non-Detection'], silent=False, verbose=True)
    #ID=[det_ID,nandet_ID]

###Implimenting Local analog calibration###
#dataset = 'O32_Grid'
out_pdf = fitspath+ '/'+dataset+'_LAC.pdf'
temp_met_ascii = fitspath+ '/'+dataset+'_temperatures_metalicity.tbl'
temp_table= asc.read(temp_met_ascii)
SN_4363 = temp_table['S/N_4363']
det_4363 = np.where((SN_4363>=3))[0]
nan_det_4363 = np.where((SN_4363<3))[0]
print det_4363
O32_all = temp_table['O32_Composite']
R23_all = temp_table['R23_Composite']
com_O_log = temp_table['com_O_log']  #This is the 12+log(OH) value
ID = temp_table['ID']
print O32_all

det_O32 = O32_all[det_4363]
det_R23 = R23_all[det_4363]
det_OH = com_O_log[det_4363]



#OH= local_analog_calibration.bian18_O32(det_O32)
#print OH
#R23_func = local_analog_calibration.bian18_R23(OH)
#xsprint R23_func

lR23 = [det_R23]
lO32 = [det_O32]
OH = [det_OH]

local_analog_calibration.main(lR23, lO32, OH, out_pdf, ID=[ID], silent=False, verbose=True)


    




###NOT IN USE###

if dataset == 'O32_Grid': outfile_single = fitspath +'single_grid_O32.npz'
    if dataset == 'R23_Grid': outfile_single = fitspath+'single_grid_R23.npz'
    
    grid_data = np.load(outfile_single)
        
    N_arr_grid = grid_data['N_arr0']
    R23_grid = grid_data['R23_grid']
    O32_grid = grid_data['O32_grid']
    T_arr = grid_data['T_arr']
        
        
    #Option to Change: Masking the night sky emission lines 
        if mask == True:
            Stack_name = 'Stacking_Masked_MasterGrid_single'+dataset+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data)
        else:
            Stack_name = 'Stacking_MasterGrid__single'+dataset+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master(fitspath, Stack_name,grid_data)
        #Outfile and pdf both use name
        print 'finished with stacking,' + Stack_name + 'pdf and fits files created'



        #Zoom_and_gauss_general
        Stack_name = Stack_name.replace('.pdf', '.fits')
        outfile_grid = fitspath + Stack_name
        print outfile_grid
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        #Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']
        binning_avg_asc = fitspath+'/'+dataset+'binning_averages.tbl'

        lineflag = np.zeros(len(wave))
        for ii in lambda0:   
            idx = np.where(np.absolute(wave - ii)<=5)[0]
            if len(idx) > 0:
                lineflag[idx] = 1
        #Option to change: Constants used as initial guesses for gaussian fit
        s=1.0
        a= 1.0
        c = 1
        
        s1=-0.3
        a1= 4.7
        s2 = 1
        a2 = -1.8

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2, binning_avg_asc)

        print 'finished gaussian fitting:,' +fitspath+'_'+dataset+'_Zoomed_Gauss_* pdfs and fits created'

        #hstack_table
        #Option to change: name of new fits file created
        intro = fitspath + 'Grid_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/Grid_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'Grid_combined_flux_table.fits'
        combine_flux_ascii = fitspath + 'Grid_combined_flux_table.tbl'
        hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_ascii)
        
        print 'Grid_combine_flux_table created'
        
        ########FIX THIS CODE##########line_ratio_plotting
        #I need to go back through and figure out what is the average and what is the composite
        line_ratio_plotting.Plotting_Data1(fitspath,dataset,combine_flux_ascii, binning_avg_asc)

       

        #R_temp_calcul
        combine_flux_ascii = fitspath + 'Grid_combined_flux_table.tbl'
        temp_m_gascii = fitspath+ '/Grid_temperatures_metalicity.tbl'
        temp_m_gfits = fitspath+ '/Grid_temperatures_metalicity.fits'
        temp_m_gpdf_name = 'Grid_Temp_Composite_Metallicity.pdf'

        R_temp_calcul.run_function(fitspath, temp_m_gascii , temp_m_gfits, temp_m_gpdf_name, combine_flux_ascii) 


        temp_table= asc.read(temp_m_gascii)
        SN_4363 = temp_table['S/N_4363']
        det_4363 = np.where((SN_4363>=3))[0]
        print det_4363'''

