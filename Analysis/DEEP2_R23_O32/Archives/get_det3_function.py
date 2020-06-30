import numpy as np
from astropy.table import Table
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
    # May limit the logR23 value further to 1.2, check the velocity dispersions of the high R23 spectra
    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (O2_ini > 0) & (O3_ini > 0) & (Hb_ini>0) & (exclude_flag==0) & (lR23_ini < 1.4))[0]  

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

    

    n2= ('ID','logR23','logO32','OII_3727_Flux_Gaussian','O3_Flux_Gaussian','HGAMMA_Flux_Gaussian',
         'HDELTA_Flux_Gaussian', 'OIII_4363_Flux_Gaussian','OIII_4958_Flux_Gaussian','OIII_5007_Flux_Gaussian',
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
        mass_valid_idx_tab = np.load('/Users/reagenleimbach/Desktop/Zcalbase_gal/From_Caroline/'
                                     'mass_bin_hbeta_revised_75_112_113_300_600_1444_1444.npz')
        mass_LHbeta_file = np.load('/Users/reagenleimbach/Desktop/Zcalbase_gal/From_Caroline/'
                                   'revised_75_112_113_300_600_1444_1444_mass_SFR_data.npz')
        
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
