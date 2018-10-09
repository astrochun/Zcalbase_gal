#This function runs the entire process start to finish
#Not completed

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
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

from Zcalbase_gal.Analysis.DEEP2_R23_O32 import Binning_and_Graphing_MasterGrid, Stackboth_MasterGrid, zoom_and_gauss_general, hstack_tables,  adaptivebinning, Stacking_voronoi, R_temp_calcul  #line_ratio_plotting,

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method_925/'
fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32]

lambda0 =[3726.16, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4363.21, 4861.32, 4958.91, 5006.84]
 
line_type = ['Oxy2', 'Balmer', 'Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']

line_name = ['OII_3727','H_9','NeIII','HeI','HEPSIL', 'HDELTA', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']   # 'NII_6548', 'HALPHA', 'NII_6584', 'SII_6717', 'SII_6730'] #'H_10' line not fitted because the 'x0 is infeasible' error occurred. In future go back and change para_bounds so that the line canbe fit  4340.46; this is the reason why we have two lists of xcoor and lambda0



def get_det3():
    for ii in range(1,5):
        file1 = fitspath_ini+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])
    #fits.writeto(fitspath_ini+'f3_0716/DEEP2_Fields_combined', data0)
        #print 'data0 : ', len(data0)
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


    data3 = data0[det3]

    O2_det3 =O2[det3]
    O3_det3 =O3[det3]
    Hb_det3 =Hb[det3]
    
    return R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3



#redue calling for look at hb


def run_R23_O32_analysis(dataset,mask='None'):
    #dataset options: Grid, Voronoi10,Voronoi14
    
    
    R23, O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3 = get_det3()
    if dataset == 'Grid':
        #Binning and Graphing MasterGrid
        #Options to Change: Bin Size
        pdf_pages1 = PdfPages(fitspath+'R23_O32_bin025_scatter_and_hexbin_MasterGrid.pdf')
        outfile1 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'
        R23_bin = 0.25
        O32_bin = 0.25
        binstr = 025
        
        Binning_and_Graphing_MasterGrid.making_Grid(fitspath, pdf_pages1,outfile1,R23,O32, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3, R23_bin, O32_bin)
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
        
        #Option to Change: Masking the night sky emission lines 
        if mask == True:
            Stack_name = 'Stacking_Wave_vs_Spect1D_Masked_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(fitspath, Stack_name,grid_data)
        else:
            Stack_name = 'Stacking_Wave_vs_Spect1D_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master(fitspath, Stack_name,grid_data)

        #Outfile and pdf both use name
        print 'finished with stacking,' + Stack_name + 'pdf and fits files created'



        #Zoom_and_gauss_general
        #outfile_grid = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/' + Stack_name
        Stack_name = Stack_name.replace('.pdf', '.fits')
        outfile_grid = fitspath + Stack_name
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']


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

        zoom_and_gauss_general.zm_general(dataset, fitspath, stack2D, header, wave, lineflag, Spect_1D, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2)#, N_gal_array=N_arr_grid, R_23_array=R23_grid, O_32_array=O32_grid)

        print 'finished gaussian fitting:,' +fitspath+'_'+dataset+'_Zoomed_Gauss_* pdfs and fits created'

        #hstack_table
        #Option to change: name of new fits file created
        intro = fitspath + 'Grid_Average_R23_O32_Values.tbl' 
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath +'/Grid_flux_gaussian_*.tbl') 
        combine_flux_table = fitspath + 'Grid_combined_flux_table.fits'
        hstack_tables.h_stack(fitspath, table_files, asc_intro, combine_flux_table)

        print 'Grid_combine_flux_table created'

        #line_ratio_plotting
        #I need to go back through and figure out what is the average and what is the composite
        #line_ratio_plotting.Plotting_Data1()


        #R_temp_calcul
        combine_flux_table = fitspath + 'Grid_combined_flux_table.fits'
        combine_fits, header_combine_fits = fits.getdata(combine_flux_table, header = True)
        print header_combine_fits
        out_ascii = fitspath+ '/Grid_temperatures_metalicity.tbl'
        out_fits = fitspath+ '/Grid_temperatures_metalicity.fits'
        pdf_name = 'Grid_Temp_Composite_Metallicity'

        R_temp_calcul.run_function(combine_fits, header_combine_fits, fitspath, out_ascii, out_fits, pdf_name)

        





'''


    if dataset == 'Voronoi10' or 'Voronoi14':
        #Adaptive Binning 
        #Options to Change: Signal/Noise Size
        if dataset == 'Voronoi10':
            sn_size = 10
            txt_file = fitspath + 'voronoi10_2d_binning_outputs.txt'
            asc_table1 = fitspath+'voronoi10_binning_averages.tbl'
            asc_table2 = fitspath+'voronoi10_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_10.tbl'
            
        if dataset == 'Voronoi14':
            sn_size = 14
            txt_file = fitspath + 'voronoi14_2d_binning_outputs.txt'
            asc_table1 = fitspath+'voronoi14_binning_averages.tbl'
            asc_table2 = fitspath+'voronoi14_2d_binning_datadet3.tbl' #used to be called fitspath+'voronoi_2d_binning_output_14.tbl'


        adaptivebinning.voronoi_binning_DEEP2(fitspath, sn_size,txt_file, asc_table1, asc_table2)


        #Stacking_voronoi
        #Option to Change: 
        if dataset == "Voronoi10": outfile = fitspath+ ''
        if dataset == "Voronoi14": outfile = fitspath + ''
        
        voronoi_data = np.load(outfile)
        

        #Option to Change: Masking the night sky emission lines 
        if mask == True:
            Stack_name = 'Stacking_Wave_vs_Spect1D_Masked_Mastergrid_voronoi'+str(binstr)+'.pdf'
            Stacking_voronoi.run_Stacking_Master_mask(fitspath, Stack_name)
        else:
            Stack_name = 'Stacking_Wave_vs_Spect1D_MasterGrid_voronoi'+str(binstr)+'.pdf'
            Stacking_voronoi.run_Stacking_Master(fitspath, Stack_name)

        #Outfile and pdf both use name
        print 'finished with stacking,' + Stack_name + 'pdf and fits files created'



        #Zoom_and_gauss_general

        outfile_grid = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/' + Stack_name
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']

        #Option to change: Constants used as initial guesses for gaussian fit
        s=1.0
        a= 1.0
        c = 1
        
        s1=-0.3
        a1= 4.7
        s2 = 1
        a2 = -1.8

        zoom_and_gauss_general.zm_general(stack2D, header, wave, Spect_1D, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2)

        print 'finished gaussian emission fitting pdf and tables created'

        #hstack_table
        #Option to change: name of new fits file created
        intro = fitspath + 'Voronoi_Average_R23_O32_Values.tbl'
        asc_intro = asc.read(intro)
        table_files = glob.glob(fitspath+ 'Voronoi_flux_gaussian_*.tbl' )
        combine_flux_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_combined_flux_table.fits'
        hstack_tables.h_stack(table_files, asc_intro, combine_flux_table)

        print 'Grid_combine_flux_table created'

        #line_ratio_plotting
        #I need to go back through and figure out what is the average and what is the composite

        #R_temp_calcul
        combine_fits = fits.getdata(combine_flux_table, header = True)
        out_ascii = fitspath+ '/Grid_temperatures_metalicity.tbl'
        out_fits = fitspath+ '/Grid_temperatures_metalicity.fits'
        pdf_name = 'Grid_Temp_Composite_Metallicity'

        R_temp_calcul.run_function(out_ascii, out_fits, pdf_name)









    if dataset == 'Voronoi14': 
         #Binning and Graphing MasterGrid
        #Options to Change: Bin Size
        pdf_pages1 = PdfPages(fitspath+'R23_O32_bin025_scatter_and_hexbin_MasterGrid.pdf')
        outfile1 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'
        R23_bin = 0.25
        O32_bin = 0.25
        binstr = 025
        
        Binning_and_Graphing_MasterGrid.making_Grid(pdf_pages1,outfile1, O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3, R23_bin, O32_bin)
        print 'made Arrays_R23O32bin025MasterGrid.npz, scatter_and_hexbin, testmastergrid(need to find if this is used anywhere)'
        print 'finished Binning_and_Graphing_MasterGrid'


        #Stackboth_MasterGrid
        #Option to Change: Bin size 
        outfile01 = fitspath+ 'Arrays_R23O32bin01MasterGrid.npz'
        outfile025 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz'
        if R23_bin == 0.25:
            grid_data = np.load(outfile025)
        else: grid_data = np.load(outfile01)

        #Option to Change: Masking the night sky emission lines 
        if mask == True:
            Stack_name = 'Stacking_Wave_vs_Spect1D_Masked_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master_mask(Stack_name)
        else:
            Stack_name = 'Stacking_Wave_vs_Spect1D_MasterGrid_bin'+str(binstr)+'.pdf'
            Stackboth_MasterGrid.run_Stacking_Master(Stack_name)

        #Outfile and pdf both use name
        print 'finished with stacking,' + Stack_name + 'pdf and fits files created'



        #Zoom_and_gauss_general

        outfile_grid = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/' + Stack_name
        stack2D, header = fits.getdata(outfile_grid, header=True)
        wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
        Spect_1D = fits.getdata(outfile_grid)
        dispersion = header['CDELT1']

        #Option to change: Constants used as initial guesses for gaussian fit
        s=1.0
        a= 1.0
        c = 1
        
        s1=-0.3
        a1= 4.7
        s2 = 1
        a2 = -1.8

        zoom_and_gauss_general.zm_general(stack2D, header, wave, Spect_1D, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2)

        print 'finished gaussian fitting,' ++ 'created'

        #hstack_table
        #Option to change: name of new fits file created
        intro = fitspath + #fill in name
        asc_into = asc.read(intro)
        table_files = glob.glob(fitspath )+  #fill in name)
        combine_flux_table = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_combined_flux_table.fits'
        hstack_tables.h_stack(table_files, asc_intro, combine_flux_table)

        print 'Grid_combine_flux_table created'

        #line_ratio_plotting
        #I need to go back through and figure out what is the average and what is the composite

        #R_temp_calcul
        combine_fits = fits.getdata(combine_flux_table, header = True)
        out_ascii = fitspath+ '/Grid_temperatures_metalicity.tbl'
        out_fits = fitspath+ '/Grid_temperatures_metalicity.fits'
        pdf_name = 'Grid_Temp_Composite_Metallicity'

        R_temp_calcul.run_function(out_ascii, out_fits, pdf_name)'''
