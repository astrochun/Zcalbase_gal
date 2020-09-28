import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from astropy.table import Table
import os
from os.path import exists, join
import glob
from datetime import date

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
