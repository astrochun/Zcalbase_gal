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

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import Binning_and_Graphing_MasterGrid, Stackboth_MasterGrid, zoom_and_gauss_general, hstack_tables,  adaptivebinning, Stacking_voronoi, R_temp_calcul, line_ratio_plotting

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Gridmethod_1130/'
fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]

lambda0 =[3726.16, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]  #, '4101.73, 4363.21, 4861.32']
        #[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]

line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single', 'Balmer', 'Balmer', 'Balmer', 'Single', 'Balmer','Single', 'Single']

line_name = ['OII_3727','H_10', 'H_9','NeIII','Hzeta','HeI','3967', 'HEPSIL', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']  

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

    
    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini>0) & (exclude_flag==0))[0]


    data3 = data0[det3]

    R23 = R23_ini[det3]
    O32 = O32_ini[det3]
    Hb  = Hb_ini[det3]
    O2  = O2_ini[det3]
    O3  = O3_ini[det3]
    SNR2 =SNR2_ini[det3]
    SNR3 =SNR3_ini[det3]
    SNRH =SNRH_ini[det3]
    
    O2_det3 =O2_ini[det3]
    O3_det3 =O3_ini[det3]
    Hb_det3 =Hb_ini[det3]
    
    return R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3



#redue calling for look at hb


def run_R23_O32_analysis(dataset,mask='None'):
    #dataset options: Grid, Voronoi10,Voronoi14, Voronoi20
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3 = get_det3()
    if dataset == 'Grid':
        #Binning and Graphing MasterGrid
        #Options to Change: Bin Size
        pdf_pages_hex = PdfPages(fitspath+'R23_O32_bin025_scatter_and_hexbin_MasterGrid.pdf')
        outfile1_npz = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'
        R23_bin = 0.25
        O32_bin = 0.25
        binstr = 025
        
        Binning_and_Graphing_MasterGrid.making_Grid(fitspath, pdf_pages_hex,outfile1_npz,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3, R23_bin, O32_bin)
        print 'made Arrays_R23O32bin025MasterGrid.npz, scatter_and_hexbin, testmastergrid(need to find if this is used anywhere)'
        print 'finished Binning_and_Graphing_MasterGrid'
        

        #Stackboth_MasterGrid
        #Option to Change: Bin size 
        outfile01 = fitspath+ 'Arrays_R23O32bin01MasterGrid.npz'
        outfile025 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz' #this file has the average R23 and O32 values for grid method
        if R23_bin == 0.25:
            grid_data = np.load(outfile025)
        else: grid_data = np.load(outfile01)
        N_arr_grid = grid_data['N_arr0']
        R23_grid = grid_data['R23_grid']
        O32_grid = grid_data['O32_grid']
        T_arr = grid_data['T_arr']
        
        
        #Option to Change: Masking the night sky emission lines 
        if mask == True:
            Stack_name = 'Stacking_Masked_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(det3, data3, fitspath,fitspath_ini, dataset, Stack_name,grid_data)
        else:
            Stack_name = 'Stacking_MasterGrid_bin'+str(binstr)+'.pdf'
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

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2,tab = binning_avg_asc)

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
        print det_4363














###Voronoi###

    if dataset == 'Voronoi10' or 'Voronoi14' or 'Voronoi20':
        #Adaptive Binning 
        #Options to Change: Signal/Noise Size
        if dataset == 'Voronoi10':
            sn_size = 10
            txt_file = fitspath + 'voronoi10_2d_binning_outputs.txt'
            asc_table1 = fitspath+'voronoi10_binning_averages.tbl' #used to be called asc_tab_fill in name
            asc_table2 = fitspath+'voronoi10_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_10.tbl'
            
        if dataset == 'Voronoi14':
            sn_size = 14
            txt_file = fitspath + 'voronoi14_2d_binning_outputs.txt'
            asc_table1 = fitspath+'voronoi14_binning_averages.tbl'
            asc_table2 = fitspath+'voronoi14_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_14.tbl'

        if dataset == 'Voronoi20':
            sn_size = 20
            txt_file = fitspath + 'voronoi20_2d_binning_outputs.txt'
            asc_table1 = fitspath+'voronoi20_binning_averages.tbl'
            asc_table2 = fitspath+'voronoi20_2d_binning_datadet3.tbl' 

        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2, R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3)

        #Stacking_voronoi
        #Option to Change: 
        if dataset == "Voronoi10": voronoi_data = fitspath+ 'voronoi10_2d_binning_datadet3.tbl'
        if dataset == "Voronoi14": voronoi_data = fitspath + 'voronoi14_2d_binning_datadet3.tbl'
        if dataset == "Voronoi20": voronoi_data = fitspath + 'voronoi20_2d_binning_datadet3.tbl'

        print '### outfile for datadet3: '+voronoi_data
        
        #Option to Change: Masking the night sky emission lines
        ######Check to make sure tables are going into right places#####
        if mask == True:
            Stack_name = 'Stacking'+dataset+'_output.pdf'
            Stacking_voronoi.run_Stacking_Master_mask(fitspath_ini, fitspath, voronoi_data, det3, asc_table1, Stack_name)
        else:
            Stack_name = 'Stacking'+dataset+'_output.pdf'
            Stacking_voronoi.run_Stacking_Master(fitspath_ini, fitspath, voronoi_data, det3, Stack_name)

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
        s=1.0
        a= 1.0
        c = 1
        
        s1=-0.3
        a1= 4.7
        s2 = 1
        a2 = -1.8
        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave,lineflag, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2,tab=asc_table1)

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

        R_temp_calcul.run_function(fitspath,temp_met_ascii,temp_met_fits, pdf_name_temp_met, combine_flux_ascii)
  

        temp_table= asc.read(temp_met_ascii)
        SN_4363 = temp_table['S/N_4363']
        det_4363 = np.where((SN_4363>=3))[0]
        print det_4363

