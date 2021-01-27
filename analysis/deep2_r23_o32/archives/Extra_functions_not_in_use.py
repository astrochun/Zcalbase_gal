from Zcalbase_gal.analysis.deep2_r23_o32 import get_det3
#run_Stacking_Master_mask()
#run_Stacking_Master()

#Plotting Zoomed in on 4363


def zoom_plot_4363():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'/Results_Graphs_Arrays_Etc/Arrays/Stacking_Wave_vs_Spect1D_Masked_MasterGrid_Average_bin025.fits')
    name = 'Stacking_Zoomed_4363.pdf'
    pdf_pages = PdfPages(fitspath+name)
    nrows = 4
    ncols = 4
    R23_grid = grid_data['R23_grid']
    O32_grid = grid_data['O32_grid']
    #index= grid_data['T_arr'][rr,rr]
    print(Spect_1D.shape)
    
    for rr in range(Spect_1D.shape[0]):

        row = rr / nrows % ncols
        col = rr % ncols
        print(row, col)
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1, x2  = 4200, 4500

        y_smooth = movingaverage_box1D(Spect_1D[rr], 5, boundary='extend')

        t_ax.plot(wave, y_smooth, linewidth=0.5)
        t_ax.set_xlim([x1,x2])
        #t_ax.set_xlabel('Wavelength')
        #t_ax.set_ylabel('Spectra 1D')
        for x in xcoor: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

        #txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 = r'R_23%.3f O32 %.3f\n' % (grid_data['R23_grid'][rr], grid_data['O32_grid'][rr])  #$\overline{x}$:$\overline{y}$:
        #txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr])
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')

        str0 = 'R23=%.1f O32=%.1f ' % (R23_grid[rr], O32_grid[rr])# N=%i  len(index))
        t_ax.annotate(str0, (0.05,0.95), xycoords='axes fraction', ha='left', va='top', weight='bold')
               

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        if col == 0:
             t_ax.set_ylabel('Spect_1D')
        else: t_ax.set_yticklabels([])
    
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.2, right=0.98, bottom=0.06, top=0.95,
                            hspace=0.0)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(11,8)
            plt.draw()
            fig.savefig(pdf_pages, format='pdf')
        
    

                                             
    #endfor
    pdf_pages.close()
    print 'Done!' 






#NOT CURRENTLY IN USE
#Plotting Zoomed in on 4363
def zoom_plot_4363():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
    name = 'Stacking_Voronoi_Zoomed_4363.pdf'
    pdf_pages = PdfPages(fitspath+name)
    nrows = 4
    ncols = 4
    
    for rr in range(Spect_1D.shape[0]):

        row = rr / nrows % ncols
        col = rr % ncols
        print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1, x2  = 4200, 4500

        y_smooth = movingaverage_box1D(Spect_1D[rr], 5, boundary='extend')

        t_ax.plot(wave, y_smooth, linewidth=0.5)
        t_ax.set_xlim([x1,x2])
        #t_ax.set_xlabel('Wavelength')
        #t_ax.set_ylabel('Spectra 1D')
        for x in xcoor: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

        txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 += 'R_23%.3f O32 %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])  #$\overline{x}$:$\overline{y}$:
        txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr])
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        if col == 0:
             t_ax.set_ylabel('Spect_1D')
        else: t_ax.set_yticklabels([])
    
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.2, right=0.98, bottom=0.06, top=0.95,
                            hspace=0.0)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(11,8)
            plt.draw()
            fig.savefig(pdf_pages, format='pdf')
        
    

                                             
    #endfor
    pdf_pages.close()
    print 'Done!'

#Plotting Zoomed in on 5007    
def zoom_plot_5007():
    image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
    name = 'Stacking_Voronoi_Zoomed_5007.pdf'
    pdf_pages = PdfPages(fitspath+name)
    nrows = 4
    ncols = 4
    
    for rr in range(Spect_1D.shape[0]):

        row = rr / nrows % ncols
        col = rr % ncols
        print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1, x2  = 4950, 5050

        y_smooth = movingaverage_box1D(Spect_1D[rr], 5, boundary='extend')

        t_ax.plot(wave, y_smooth, linewidth=0.5)
        t_ax.set_xlim([x1,x2])
        #t_ax.set_xlabel('Wavelength')
        #t_ax.set_ylabel('Spectra 1D')
        for x in xcoor: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

        txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
        txt0 += 'R_23%.3f O32 %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])  #$\overline{x}$:$\overline{y}$:
        txt0 += 'S/N: %.3f  Scale: %.3f N: %.3f\n' % (asc_tab['sn'][rr], asc_tab['scale'][rr], asc_tab['area'][rr])
       
        t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '6')

        if row != nrows-1:
            t_ax.set_xticklabels([]) 
        else: t_ax.set_xlabel('Wavelength')

        if col == 0:
             t_ax.set_ylabel('Spect_1D')
        else: t_ax.set_yticklabels([])
    
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.2, right=0.98, bottom=0.06, top=0.95,
                            hspace=0.0)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(11,8)
            plt.draw()
            fig.savefig(pdf_pages, format='pdf')
        
    

                                             
    #endfor
    pdf_pages.close()
    print 'Done!'



'''#Initialize Error Propogation
    flux_g_err = np.zeros(Spect_1D.shape[0])
    flux_s_err = np.zeros(Spect_1D.shape[0])
    sigma_err = np.zeros(Spect_1D.shape[0])
    median_err = np.zeros(Spect_1D.shape[0])
    norm_err = np.zeros(Spect_1D.shape[0])
    RMS_err = np.zeros(Spect_1D.shape[0])
    SN_err = np.zeros(Spect_1D.shape[0])'''


            '''#Error Propogation 
            flux_g_err[rr] = error_prop_chuncodes(flux_g,1)
            flux_s_err[rr] = error_prop_chuncodes(flux_s,1)
            sigma_err[rr] =  error_prop_chuncodes(o1[1],1)
            median_err[rr] = error_prop_chuncodes(o1[3],1)
            norm_err[rr] =   error_prop_chuncodes(max0,1)
            RMS_err[rr] =    error_prop_chuncodes(ini_sig1,1)
            SN_err[rr] =     error_prop_chuncodes(flux_s/ini_sig1,1)'''

            '''if row != nrows-1:
                t_ax.set_xticklabels([]) 
            else: t_ax.set_xlabel('Wavelength')

            #if row == nrows-1: t_ax.set_xlabel('Wavelength')

            if col == 0:
                t_ax.set_ylabel('Spect_1D')
                t_ax.set_yticklabels([])  #sets y-tick labels 
    
            if rr == Spect_1D.shape[0]-1 and rr % (nrows*ncols) != nrows*ncols-1:
                for jj in range(col+1, ncols): ax_arr[row,jj].axis('off')
            for kk in range(row+1, nrows):
                for zz in range(ncols): ax_arr[kk,zz].axis('off')

            xlabels = t_ax.get_xticklabels()
            #print 'xlabels : ', xlabels
            for yy in range(rr,rr-ncols,-1):
                y_col = yy % ncols
                y_row = yy / nrows % ncols
                ax_arr[y_row,y_col].set_xlabel('Wavelength')
                ax_arr[y_row,y_col].set_xticklabels(xlabels)'''


'''
    #Error Propogation 
    flux_g_err = error_prop_chuncodes(flux_g_array,1)
    flux_s_err = error_prop_chuncodes(flux_s_array,1)
    sigma_err =  error_prop_chuncodes(sigma_array,1)
    median_err = error_prop_chuncodes(median_array,1)
    norm_err =   error_prop_chuncodes(norm_array,1)
    RMS_err =    error_prop_chuncodes(RMS_array,1)
    SN_err =     error_prop_chuncodes(SN_array,1)'''


'''
        print 'Starting fits'
        #header['Comment'] = 'This fits file contains the calculated and observed fluxs for emission line stated in file name'
        c1 = fits.Column(name='Flux_Gaussian'+line_name, format= 'K', array=flux_g_array)
        c2 = fits.Column(name='Flux_Observed'+line_name, format= 'K', array=flux_s_array)
        c3 = fits.Column(name='Sigma'+line_name,  format= 'K', array=sigma_array)
        c4 = fits.Column(name='Median'+line_name,  format= 'K', array=median_array)
        c5 = fits.Column(name='Norm'+line_name, format= 'K',  array=norm_array)
        c6 = fits.Column(name='RMS'+line_name, format= 'K',  array=ini_sig1)

        #fits_cols = fits.ColDefs([c1,c2,c3,c4,c5])
        f_tab= fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])
        
        f_tab.writeto(fitspath+'/Mar19_outputs/Flux_Outputs'+line_name+'.fits', overwrite= True)'''


#Unnecessary Codes:
 #Calculating Gauss 
 #if keyword == 'Single':
 #    gauss1 = gauss(x0,*o1)
 #if keyword == 'Balmer':
 #    gauss1 = double_gauss(x0,*o1)
 #if keyword == 'Oxy2':
 #    gauss1 = oxy2_gauss(x0,*o1)



 #np.where((x0-workingwave)<=100 & ((line_flag=0))
 #+-100
 #np.zero(len(x0))

 #lyetal(2014)


 

'''


#Fits Files 
#fit_table_files = glob.glob('/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar2_outputs/Flux_Outputs*.tbl')
outfile = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar2_outputs/combined_flux_table.fits'
hstacking.write(outfile, format='fits', overwrite=True)

#fits.writeto('/Users/reagenleimbach/Desktop/Zcalbase_gal/Fits_File_2_23/combined_flux_table.fits', hstacking,
#             overwrite=True)


'''


'''
#Voronoi10
intro = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi_asc_table_Average_R23_O32_Values.tbl'
asc_intro = asc.read(intro)
table_files= glob.glob('/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi_asc_table_flux_gaussian_*.tbl')

#Grid
intro = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_asc_table_Average_R23_O32_Values.tbl'
asc_intro = asc.read(intro)
table_files= glob.glob('/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_asc_table_flux_gaussian_*.tbl')

#Voronoi14
intro = '/Users/reagenleimbach/Desktop/Zcalbase_gal/MAY25_Voronoi14/Voronoi14_asc_table_Average_R23_O32_Values.tbl'
asc_intro = asc.read(intro)
table_files= glob.glob('/Users/reagenleimbach/Desktop/Zcalbase_gal/MAY25_Voronoi14/Voronoi14_asc_table_flux_gaussian_*.tbl')'''


#fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_104/'



#single_vs_averaged = fitspath +'single_vs_averaged_emissions.pdf'
#A1_attempt2 = fitspath +'A1_attempt2.pdf'
#A2_attempt2 = fitspath +'A2_attempt2.pdf'
#comparing_composite_values = fitspath +'comparing_composite_values.pdf'
'''
#Raw R23 and O32 (xnode is R23 and ynode is O32) these are average values calculated by voronoi code 10
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi.tbl'
data0 = asc.read(raw)
R23_raw_voronoi = data0['xBar']
O32_raw_voronoi = data0['yBar']'''
'''
#Raw R23 and O32 (xnode is R23 and ynode is O32) these are average values calculated by voronoi code 14
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_14.tbl'
data0 = asc.read(raw)
R23_raw_voronoi = data0['xBar']
O32_raw_voronoi = data0['yBar']

#Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Voronoi Outputs10: R_23 and O_32 values for all galaxies with a column specifying bins 
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar21_voronoi_2d_binning_output.txt'
voronoi = np.genfromtxt(outfilevoronoi)

#Voronoi Outputs14: R_23 and O_32 values for all galaxies with a column specifying bins 
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output_14.txt'
voronoi = np.genfromtxt(outfilevoronoi)

#Grid Inputs and Outputs
raw = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_asc_table_Average_R23_O32_Values.tbl' 
data0 = asc.read(raw)
R23_raw_voronoi = data0['R_23_Average']
O32_raw_voronoi = data0['O_32_Average']

spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/grid_combined_flux_table.tbl'
data1 = asc.read(spectral)'''




####The rest of this code is not fully developed###
#RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'

def average_vals_bin2(): 
    '''for ii in range(1,5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

        #print 'data0 : ', len(data0)
    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']'''

    O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3, O2_det3, O3_det3, Hb_det3 = get_det3()
    #O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3= get_det3()

    O2_det3=O2[det3]
    O3_det3=O3[det3]
    Hb_det3=Hb[det3]
    
    #Calculating R23 and O32 for the orginal data
    '''R23_original = np.zero(len(O2_det3))
    O32_original = np.zero(len(O2_det3))
    for ii in len(O2_det2): 
        R23 = (O2_det3[ii]+O3_det3[ii])/Hb_det3[ii]
        O32 = O3_det3[ii]/O2_det3[ii]
        R23_original[ii] = R23
        O32_original[ii] = O32'''
       
    #Composite Spectrum
    lR23 = voronoi[:,2]
    lO32 = voronoi[:,3]
    binnum = voronoi[:,4]
    binnum0= list(set(binnum))
    
    R23_original = []  #np.zeros(len(binnum0), dtype=np.float64)
    O32_original = []  #np.zeros(len(binnum0), dtype=np.float64)
    #print type(avg_l_R23), type(avg_l_O32)
    #print avg_l_R23
    for ii in binnum0:
        idx = np.array([xx for xx in range(len(det3)) if binnum[xx] == ii])
            
        O2_original = np.mean(O2_det3[idx])/1e-17
        O3_original = np.mean(O3_det3[idx])/1e-17
        Hb_original = np.mean(Hb_det3[idx])/1e-17
        #print 'Avg_O3', avg_O3
        list_R23 = np.log10((O2_original +O3_original)/Hb_original)
        list_O32 = np.log10(O3_original/O2_original)
        #print list_R23, type(list_R23), list_O32, type(list_O32)

        R23_original.append(list_R23)
        O32_original.append(list_O32)
    Plotting_Data1(R23_original, O32_original)
 






#The rest of the in this file does different parts of the above two functions. The functions above are the products of what is below. 

def A1_attempt_two_line_plotting():
    #When I loop over the for loop, the data from all previous loops also graphs, how do I fit that? 
    com_R23 = data1['R_23_Average']
    com_O32 = data1['O_32_Average']
    
    binNum = voronoi[:,4]
    
    
    binNum0 = list(set(binNum))
    R23_arr = np.zeros(len(binNum0))
    O32_arr = np.zeros(len(binNum0))

    pdf_pages= PdfPages(A1_attempt2)
    for ii in range(len(binNum0)):
        fig, ax_arr = plt.subplots()
        random_bin = binNum0[ii]
        bin_idx = [xx for xx in range(len(voronoi)) if binNum[xx] == random_bin]

        R23_arr= voronoi[:,1][bin_idx]
        O32_arr= voronoi[:,2][bin_idx]

        c_R23 = com_R23[ii]
        c_O32 = com_O32[ii]
     
        ax_arr.scatter(R23_arr, O32_arr)
        ax_arr.plot(c_R23, c_O32, '*',markersize=12)
        ax_arr.set_xlabel('R23')
        ax_arr.set_ylabel('O32')
        fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()
    
def A2_attempt_two_line_plotting():
    binnum = list(range(0,27))
    pdf_pages= PdfPages(A2_attempt2)
    scalefact = 1e-17
#Stuff for every galaxies
    for ii in range(1,5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

        #print 'data0 : ', len(data0)
    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']
     
    SNR2 = data0['OII_SNR']
    SNR3 = data0['OIIIR_SNR']
    SNRH = data0['HB_SNR']
    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                    (O2 > 0) & (O3 > 0) & (Hb>0))[0]


    data3 = data0[det3]
    O2 = data3['OII_FLUX_MOD']/scalefact
    O3 = 1.33*data3['OIIIR_FLUX_MOD']/scalefact
    Hb = data3['HB_FLUX_MOD']/scalefact

   
    com_O2 = data1['OII_3727_Flux_Observed']
    com_O3 = data1['OIII_4958_Flux_Observed']
    com_Hb = data1['HBETA_Flux_Observed']

    

    '''for ii in binnum:
        fig, ax_arr = plt.subplots()
        ax_arr.scatter(O2,O3) #These numbers are not binned and so I need to figure out how to call the binning 
        ax_arr.plot(com_O2[ii],com_O3[ii], '*', markersize=12)
        fig.savefig(pdf_pages, format='pdf')

    for ii in range(binnum):
        

    pdf_pages.close()'''




def deal_with_later():
#compare the averages for the individual lines that were calculated on Friday with the large table from zm_general
    O2_observed = np.log10(data1['OII_3727_Flux_Observed'])
    O3_observed = np.log10(data1['OIII_5007_Flux_Observed'])
    Hb_observed = np.log10(data1['HBETA_Flux_Observed'])
    #print O2_observed
    binnum = data1['N_Galaxies']
    pdf_pages= PdfPages(single_vs_averaged)
    print avg_O2, avg_O2.size, avg_O3, avg_O3.size

    #O2
    '''row = rr / nrows % ncols
        col = rr % ncols
        #print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]'''
    nrows = 4
    ncols = 4 
    
    fig, ax_arr = plt.subplots()
    
    for ii in range(O2_observed.shape[0]): 
   
      ax_arr.scatter(O2_observed, O3_observed, label= 'O2_observed')
      ax_arr.plot(avg_O2, avg_O3)
      #ax_arr.plot(aver_O2, aver_O3)
      #ax_arr.legend(loc=0)
      #for bb in range(len(data1)):
      #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
      #ax_arr.set_xlabel()
      ax_arr.set_ylabel('O2_Flux')

      #ax_arr.plot([0.0,1.3], [0.0,1.3], 'k-')

      fig.savefig(pdf_pages, format='pdf')

    '''

    #O3
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(binnum,O3_observed, label= 'O3')
    ax_arr.scatter(binnum,avg_O3)
    ax_arr.legend(loc=0)
    #for bb in range(len(data1)):
        #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
    ax_arr.set_xlabel('N_Galaxies')
    ax_arr.set_ylabel('O3_Flux')
    fig.savefig(pdf_pages, format='pdf')

    #Hb
    fig, ax_arr = plt.subplots()
    ax_arr.scatter(binnum,Hb_observed, label= 'Hb')
    ax_arr.scatter(binnum,avg_Hb)
    ax_arr.legend(loc=0)
    #for bb in range(len(data1)):
        #ax_arr.annotate(str(binnum[bb]), (R23_raw[bb],R23_composite[bb]), xycoords='data')
    ax_arr.set_xlabel('N_Galaxies')
    ax_arr.set_ylabel('Hb_Flux')
    fig.savefig(pdf_pages, format='pdf')
'''
    pdf_pages.close()







#R23 = ([OII] + [OIII]4959 +[OIII]5007)/H_BETAemission
#O23 = ([OIII]4959 +[OIII]5007)/[OII]
#First you need to find the R23 and O32 values for every line for all emission plots (1-27) using the ascii table. Use the spectral fluxs you calculated
#Then calculate the R23 and O32 values from the raw data(x_node and y_node values yes from Stacking_Voronio) 
#Finally generate plots with the raw R23 and the spectral R23 (same for O32) == you should find a one to one ratio: plot the log vs. log 



#Random Bin xBar and yBar Calculations
#Calculate the average R23 and O32 values for a random bin == import voronoi_2d_binning_output.txt file and for a random bin fidn the number of galaxies and sum??

def checking_R23_O32():
    outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
    voronoi = asc.read(outfilevoronoi)
    binNum = voronoi['col3']

    binNum0 = list(set(binNum))
    R23_arr = np.zeros(len(binNum0))
    O32_arr = np.zeros(len(binNum0))
    
    for ii in range(len(binNum0)):
        random_bin = binNum0[ii]
        bin_idx = [xx for xx in range(len(voronoi)) if binNum[xx] == random_bin]


        R23_arr[ii] = np.average(voronoi['col1'][bin_idx])
        O32_arr[ii] = np.average(voronoi['col2'][bin_idx])
    n= ('R_23_Average', 'O_32_Average', 'bin_number')
    tab1 = Table([R23_arr, O32_arr,binNum0], names=n)
    asc.write(tab1, fitspath+'avg_calcul_values.tbl', format='fixed_width_two_line')


#initialize arrays and populate 

def comparing_tables():
    outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
    voronoi = asc.read(outfilevoronoi)
    avg_calcul_tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/avg_calcul_values.tbl'
    calcul_tab = asc.read(avg_calcul_tab)
    #data1 is the other table
    pdf_pages = PdfPages('/Users/reagenleimbach/Desktop/Zcalbase_gal/comparing_tables.pdf')
    
    binNum = voronoi['col3']
    print len(binNum)
    
    #for ii in range(len(binNum)):
    #computer calculated - checking_R23_O32 calucation 
    diff_R23 = data1['R_23_Average']- calcul_tab['R_23_Average']
    diff_O32 = data1['O_32_Average']- calcul_tab['O_32_Average']
    #if diff_R23 != 0: print diff_R23
    #if diff_O32 != 0: print diff_O32

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(diff_R23,diff_O32)
    #ax_arr.scatter(binNum[ii],diff_O32)
    #ax_arr.set_xlabel(r'Raw log($R_{23}$)')
    #ax_arr.set_ylabel(r'Observed log($R_{23}$)')
    
    '''
    fig2, ax_arr2 = plt.subplots()
    ax_arr2.scatter(,O32_observed, label= 'O32 Ratio')       
    ax_arr2.set_xlabel(r'Raw log($O_{32}$')
    ax_arr2.set_ylabel(r'Observed log($O_{32}$)')'''


    plt.draw()
    fig.savefig(pdf_pages, format='pdf')
    #print 'end for'
    pdf_pages.close()
    print 'Finished'


def save_for_later():
    lR23 = voronoi[:,2]
    lO32 = voronoi[:,3]
    binnum = voronoi[:,4]
    binnum0= list(set(binnum))
    
    avg_l_R23 = []  #np.zeros(len(binnum0), dtype=np.float64)
    avg_l_O32 = []  #np.zeros(len(binnum0), dtype=np.float64)
    avg_l_O2 = []
    avg_l_O3 = []
    avg_l_Hb = []
    #print type(avg_l_R23), type(avg_l_O32)
    #print avg_l_R23
    for ii in binnum0:
        idx = np.array([xx for xx in range(len(det3)) if binnum[xx] == ii])
            
        avg_O2 = np.mean(O2[idx])/1e-17
        avg_O3 = np.mean(O3[idx])/1e-17
        avg_Hb = np.mean(Hb[idx])/1e-17
        print 'Avg_O3', avg_O3
        list_R23 = np.log10((avg_O2 +avg_O3)/avg_Hb)
        list_O32 = np.log10(avg_O3/avg_O2)
        #print list_R23, type(list_R23), list_O32, type(list_O32)

        avg_l_R23.append(list_R23)
        avg_l_O32.append(list_O32)
        avg_l_O2.append(avg_O2)
        avg_l_O3.append(avg_O3)
        avg_l_Hb.append(avg_Hb)


def comparing_composite_values():
    comparing_composite_values = fitspath +'comparing_composite_values.pdf'
    com_R23 = data1['R_23_Average'].data
    com_O32 = data1['O_32_Average'].data

    OII = data1['OII_3727_Flux_Observed'].data
    OIII4959 = data1['OIII_4958_Flux_Observed'].data
    OIII5007 = data1['OIII_5007_Flux_Observed'].data
    H_BETA = data1['HBETA_Flux_Observed'].data
    binnum = data1['N_Galaxies'].data
    #print binnum

    pdf_pages= PdfPages(comparing_composite_values)
    R23_composite = np.log10((OII + 1.33 *OIII5007)/H_BETA)
    O32_composite = np.log10((1.33*OIII5007)/OII)

    fig, ax_arr = plt.subplots()
    for ii in range(len(binnum)):
        ax_arr.plot(R23_composite[ii], O32_composite[ii], '*')
        print 'calculated', R23_composite[ii], O32_composite[ii]
        ax_arr.plot(com_R23[ii], com_O32[ii], '+')
        print 'computer', com_R23[ii], com_O32[ii]
        ax_arr.set_xlabel('R23')
        ax_arr.set_ylabel('O32')
        fig.savefig(pdf_pages, format='pdf')
    
    pdf_pages.close()

'''
#Voronoi10
#Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi10/Voronoi_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Grid
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Grid_method/Grid_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Voronoi14
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_combined_flux_table.tbl'
data1 = asc.read(spectral)

#Ascii Table Calls 
OIII5007 = data1['OIII_5007_Flux_Observed'].data
OIII4959 = data1['OIII_4958_Flux_Observed'].data
OIII4363 = data1['OIII_4363_Flux_Observed'].data
HBETA    = data1['HBETA_Flux_Observed'].data
OII3727  = data1['OII_3727_Flux_Observed'].data
R23_avg      = data1['R_23_Average'].data
O32_avg      = data1['O_32_Average'].data
N_Galaxy = data1['N_Galaxies'].data

SN_5007       = data1['OIII_5007_S/N'].data
SN_4959       = data1['OIII_4958_S/N'].data
SN_4363       = data1['OIII_4363_S/N'].data
SN_HBETA      = data1['HBETA_S/N'].data
SN_3727       = data1['OII_3727_S/N'].data


#Fits Table Calls (This is incorrect)
OIII5007 = header['OIII_5007_Flux_Observed']
OIII4959 =  header['OIII_4958_Flux_Observed']
OIII4363 =  header['OIII_4363_Flux_Observed']
HBETA    =  header['HBETA_Flux_Observed']
OII3727  =  header['OII_3727_Flux_Observed']
R23_avg      =  header['R_23_Average']
O32_avg      =  header['O_32_Average']
N_Galaxy =  header['N_Galaxies']

SN_5007       =  header['OIII_5007_S/N']
SN_4959       =  header['OIII_4958_S/N']
SN_4363       =  header['OIII_4363_S/N']
SN_HBETA      =  header['HBETA_S/N']
SN_3727       =  header['OII_3727_S/N']

R23_composite = np.log10((OII3727 + (1.33*OIII5007))/HBETA)
O32_composite = np.log10((1.33*OIII5007)/OII3727)
'''



    '''R_value = np.zeros(len(OIII5007))
    T_e = np.zeros(len(OIII5007))
    O_s_ion = np.zeros(len(OIII5007))
    O_d_ion = np.zeros(len(OIII5007))
    com_O_log = np.zeros(len(OIII5007))
    
    for ii in range(len(OIII5007)):
        print SN_4363[ii]
        if SN_4363[ii] >= 3:
            print 'regular'
            R_value_cal= R_calculation(OIII4363[ii], OIII5007[ii], OIII4959[ii])   #, SN_4636, SN_5007, SN_495)
            T_e_cal= temp_calculation(R_value_cal)  #, R_std)
            O_s_ion_cal, O_d_ion_cal, com_O_log_cal, log_O_s_cal, log_O_d_cal = metalicity_calculation(T_e_cal,OIII5007[ii], OIII4959[ii], OIII4363[ii], HBETA[ii], OII3727[ii])

            R_value[ii] = R_value_cal
            T_e[ii] = T_e_cal
            O_s_ion[ii] = O_s_ion_cal
            O_d_ion[ii] = O_d_ion_cal
            com_O_log[ii] = com_O_log_cal
        else:
            print 'upper limit'
            R_value_up_cal= R_calculation(up_limit[ii], OIII5007[ii], OIII4959[ii])   #, SN_4636, SN_5007, SN_495)
            T_e_up_cal= temp_calculation(R_value_up_cal)  #, R_std)
            O_s_ion_up_cal, O_d_ion_up_cal, com_O_log_up_cal, log_O_s_up_cal, log_O_d_up_cal = metalicity_calculation(T_e_up_cal,OIII5007[ii], OIII4959[ii], up_limit[ii], HBETA[ii], OII3727[ii])

            R_value[ii] = R_value_up_cal
            T_e[ii] = T_e_up_cal
            O_s_ion[ii] = O_s_ion_up_cal
            O_d_ion[ii] = O_d_ion_up_cal
            com_O_log[ii] = com_O_log_up_cal'''

    
    #fill in the nan detections and refill anything less than 3 sigma in with these values
    #set it so that function requires you to enter a sigma value into the original function which I think it does but you just need to pass it all the way through

    #Ascii Table
    #out_ascii = fitspath+ '/Grid_temperatures_metalicity_asc_table.tbl'
    #out_fits = fitspath+ '/Grid_temperatures_metalicity_asc_table.fits'

B_com_R23, BR23, B_com_O32, BO32 = R23_O32_relations_BIAN()

#color_arr = [b,g,y]
    #for i in range(len(R23_composite)):
        #color_arr[i]= cmap[i]
    #print 'color_arr', color_arr
    #print len(color_arr)
    #print color_len




def R23_O32_relations_BIAN():
    x =np.zeros(100)
    y =np.zeros(100)
    for a in range(5,11,100):
        c = 138.0430 - 54.8284*a + 7.2954*(a^2) -0.32293*(a^3)
        x[a] = a #np.append(x,a)
        y[a] = c #np.append(y,c)
    t = np.zeros(100)
    z = np.zeros(100)
    for b in range(-1,1,100):
        d = 8.54 - 0.59*b
        t[b] = b #np.append(t,b)
        z[b] = d #np.append(z,d)
    #print x, y, t, z

    return x, y, t, z 

'''We will adopt a S/N of 3.  So if you get S/N < 3, you should automatically fix the 4363 flux to the 3-sigma limit.  We can use H-gamma to get the line width to determine how many pixels.'''
#log(O/H)
#error propagation
'''we have a flux and an sigma flux = flux/signa to noise
propagation distrubution function... monticarlo'''

#3D plots


#Plot O/H values and try the linear plots on the voronoi 14 values 


def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)




def R23_O32_relations_BIAN():
    x =np.zeros(100)
    y =np.zeros(100)
    for a in range(5,11,100):
        c = 138.0430 - 54.8284*a + 7.2954*(a^2) -0.32293*(a^3)
        x[a] = a #np.append(x,a)
        y[a] = c #np.append(y,c)
    t = np.zeros(100)
    z = np.zeros(100)
    for b in range(-1,1,100):
        d = 8.54 - 0.59*b
        t[b] = b #np.append(t,b)
        z[b] = d #np.append(z,d)
    #print x, y, t, z

    return x, y, t, z 

'''We will adopt a S/N of 3.  So if you get S/N < 3, you should automatically fix the 4363 flux to the 3-sigma limit.  We can use H-gamma to get the line width to determine how many pixels.'''
#log(O/H)
#error propagation
'''we have a flux and an sigma flux = flux/signa to noise
propagation distrubution function... monticarlo'''

#3D plots


#Plot O/H values and try the linear plots on the voronoi 14 values 


def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)


    '''#Histogram
    name = 'Temperature_histogram.pdf'
    pdf_pages = PdfPages(fitspath+name)
    plt.hist(valid_T_e, bins =8)
    plt.xlabel('Temperature (K)')
    #plt.set_ylabel('Spectra')
    plt.title('Preliminary Temperatures')
    pdf_pages.savefig()'''

