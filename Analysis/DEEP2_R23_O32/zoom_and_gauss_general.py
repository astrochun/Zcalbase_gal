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

#Imports Error propagation codes from chun_codes
from chun_codes import random_pdf, compute_onesig_pdf

#fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'
'''
#Grid File
outfile_grid = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Stacking_Wave_vs_Spect1D_Masked_MasterGrid_Average_bin025.fits'
stack2D, header = fits.getdata(outfile_grid, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
Spect_1D = fits.getdata(outfile_grid)
dispersion = header['CDELT1']




#Voronoi10
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output.txt'
voronoi = np.loadtxt(outfilevoronoi)
voronoi = voronoi.transpose()

stacking_vor= r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Stacking_Voronoi_masked_output.fits'

tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi.tbl'
asc_tab = asc.read(tab)

stack2D, header = fits.getdata(stacking_vor, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
dispersion = header['CDELT1']

#Voronoi14
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/voronoi_2d_binning_output_14.txt'
voronoi = np.loadtxt(outfilevoronoi)
voronoi = voronoi.transpose()

stacking_vor= r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Stacking_Voronoi_14_masked_output.fits'

tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi_14.tbl'
asc_tab = asc.read(tab)

stack2D, header = fits.getdata(stacking_vor, header=True)
wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_14_masked_output.fits')
dispersion = header['CDELT1']'''


#If 'x0' is infeasible error occurs, check the para_bound values to make sure the expected values are within the range set up upper and lower limits. 

RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'
#eventually I need to incoorparate these values into the correct function but right now they are global variables 
s=1.0
a= 1.0
c = 2.0
s1= 1.3
a1= 1.5
s2 = 1
a2 = 1.8


'''[3835.38 'SINGLE' , 3970.07????, ] #[5006.84]
#[4363.21,4958.91,5006.84, 3835.38, 3970.07, 4101.73, 4340.46, 4861.32] #[5006.84]
#['OIII_4363','OIII_4959','OIII_5007','NeIII','','HDELTA','HGAMMA']
#6548.10, 6562.80, 6583.60, 6717.42, 6730.78]'''

#lambda0 =[3727.00, 3797.90, 3835.38, 3868.74, 3888.65, 3970.07, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
 
#line_type = ['Oxy2','Balmer', 'Balmer', 'Single', 'Single', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']

#line_name = ['OII_3727','H_10','H_9','NeIII','HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007'] #, 'NII_6548', 'HALPHA', 'NII_6584', 'SII_6717', 'SII_6730']

'''
lineflag = np.zeros(len(wave))
for ii in lambda0:   
    idx = np.where(np.absolute(wave - ii)<=5)[0]
    if len(idx) > 0:
        lineflag[idx] = 1'''

def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth

def gauss(x,xbar,s,a,c):
   
    return  a*np.exp(-(x-xbar)**2/(2*s**2)) + c 

def double_gauss(x, xbar, s1, a1, c, s2, a2): 

    return a1*np.exp(-(x-xbar)**2/(2*s1**2)) + c + a2*np.exp(-(x-xbar)**2/(2*s2**2))

def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    con1 = 3728.91/3726.16
    return a1*np.exp(-(x-xbar)**2/(2*s1**2)) + c + a2*np.exp(-(x-(xbar*con1))**2/(2*s2**2)) 

def get_gaussian_fit(working_wave,x0, y0, y_norm, x_idx,line_type):

    med0 = np.median(y_norm[x_idx])
    max0 = np.max(y_norm[x_idx]) - med0
    #do we need to calculate the xbar,s2 values in the Balmer fits 
    #print med0, max0

    fail = 0
    if line_type == 'Single': 
        p0 = [working_wave, 1.0, max0, med0] #must have some reasonable values
       
        para_bound = ((working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0)),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0)))

        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)
        except ValueError:
            print 'fail'
            fail = 1
            
    if line_type == 'Balmer': #if guess > med : guess = median
        p0 = [working_wave, 1.0, max0, med0, s2, -0.5*max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -max0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0,0)

        #print 'get_gaussian_fit, p0 : ', p0

        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)
        except ValueError:
            print 'fail'
            fail = 1

    if line_type == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75*max0, med0, 1.0, max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, 0.0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0, 100.0)

        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0=p0,sigma=None, bounds = para_bound) #verbose= True)
        except ValueError:
            print 'fail'
            fail = 1
   
    if not fail:
        return o1, med0, max0
    else:
        return None, med0, max0

#calculating rms

def rms_func(wave,dispersion,lineflag, lambda_in, y0, sigma_array, scalefact):
    x_idx = np.where((np.abs(wave-lambda_in)<=100) & (lineflag==0))[0]
    #ini_sig = np.zeros(Spect_1D.shape[0])
    #for rr in range(Spect_1D.shape[0]):
    #y0 = stack2D[rr]
    sigma = np.std(y0[x_idx])
    pix =  5* sigma_array / dispersion 
    s_pix = np.sqrt(pix)

    ini_sig= s_pix * sigma * dispersion
    RMS_pix = sigma*dispersion/scalefact
    
    #print '#rms func : ', len(x_idx), pix, ini_sig
       
    '''pix =  5* sigma_array[rr] / dispersion 
    s_pix = np.sqrt(pix)

    ini_sig[rr]= s_pix * sigma * dispersion'''
   
    return ini_sig/scalefact, RMS_pix

def error_prop_chuncodes(x_arr, dx):
    #x_arr = np.random.random_integers(10,10)
    x_pdf = random_pdf(x_arr,0.5)
    return x_pdf
    
#for each individual stack
#electron temperature and the R23 and O32 values 
#Plotting Zoomed 
def zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, Spect_1D, working_wave, lambda0, lineflag,  line_type = '',outpdf='', line_name=''):
    '''image2DM, header = fits.getdata(RestframeMaster, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    Spect_1D = fits.getdata(fitspath+'Stacking_Voronoi_masked_output.fits')
    if outpdf == '':
        name = 'Stacking_Voronoi_Zoomed_Gauss_generalexperiment.pdf'
        outpdf = fitspath + name'''
    print 'tab='+tab
    asc_tab = asc.read(tab)
    pdf_pages = PdfPages(outpdf)
    nrows = 4
    ncols = 4
    x_idx = np.where((wave>=(working_wave-100)) & (wave<=(working_wave+100)))[0] 
    x0 = wave#[x_idx]
    scalefact = 1e-17

    '''narr_gr = np.asarray(N_arr_grid)
    oarr_gr = np.asarray(O32_grid)
    rarr_gr = np.asarray(R23_grid)
    print narr_gr, type(narr_gr)
    print rarr_gr, type(rarr_gr)
    print oarr_gr, type(oarr_gr)'''

    #Initializing Arrays
    if dataset == 'Grid':
        flux_g_array = np.zeros(Spect_1D.shape[0])
        flux_s_array = np.zeros(Spect_1D.shape[0])
        sigma_array = np.zeros(Spect_1D.shape[0])
        median_array = np.zeros(Spect_1D.shape[0])
        norm_array = np.zeros(Spect_1D.shape[0])
        RMS_array = np.zeros(Spect_1D.shape[0])
        SN_array = np.zeros(Spect_1D.shape[0])
        N_gal_array = np.zeros(Spect_1D.shape[0])
        R_23_array = np.zeros(Spect_1D.shape[0])
        O_32_array = np.zeros(Spect_1D.shape[0])

        
        '''outfile025 = fitspath + 'Arrays_R23O32bin025MasterGrid.npz' #this file has the average R23 and O32 values for grid method
        grid_data = np.load(outfile025)
       
        N_arr_grid = grid_data['N_arr0']
        R23_grid = grid_data['R23_grid']
        O32_grid = grid_data['O32_grid']
        #Attempt 1
        for rr in range(len(R23_grid)):
            for oo in range(len(O32_grid)):
                R_23_array = R23_grid[rr]
                O_32_array = O32_grid[oo]
                index= grid_data['T_arr'][rr,oo]
                N_gal_array = len(index)
                print R_23_array, O_32_array, N_gal_array  #3 dimensional array

        #Attempt 2
        data = np.zeros(110,3)
        for rr in range(len(R23_grid)):
            for oo in range(len(O32_grid)): 
                data0 = np.array([R23_grid[rr], O32_grid[oo], N_arr_grid[rr,oo]])
                
                print data'''
    else: 
        flux_g_array = np.zeros(Spect_1D.shape[0])
        flux_s_array = np.zeros(Spect_1D.shape[0])
        sigma_array = np.zeros(Spect_1D.shape[0])
        median_array = np.zeros(Spect_1D.shape[0])
        norm_array = np.zeros(Spect_1D.shape[0])
        RMS_array = np.zeros(Spect_1D.shape[0])
        SN_array = np.zeros(Spect_1D.shape[0])
        N_gal_array = np.zeros(Spect_1D.shape[0])
        R_23_array = np.zeros(Spect_1D.shape[0])
        O_32_array = np.zeros(Spect_1D.shape[0])

    '''#Initialize Error Propogation
    flux_g_err = np.zeros(Spect_1D.shape[0])
    flux_s_err = np.zeros(Spect_1D.shape[0])
    sigma_err = np.zeros(Spect_1D.shape[0])
    median_err = np.zeros(Spect_1D.shape[0])
    norm_err = np.zeros(Spect_1D.shape[0])
    RMS_err = np.zeros(Spect_1D.shape[0])
    SN_err = np.zeros(Spect_1D.shape[0])'''

    for rr in range(Spect_1D.shape[0]):
        #print range(Spect_1D.shape[0])
        y0 = stack2D[rr]
        y_norm = y0/scalefact

        row = rr / nrows % ncols
        col = rr % ncols
        #print row, col
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)
                #endif
       
        t_ax = ax_arr[row,col]
        
        x1 = working_wave-100
        x2  = working_wave+100

        y_smooth = movingaverage_box1D(Spect_1D[rr]/scalefact, 2, boundary='extend')
      

        o1, med0, max0  = get_gaussian_fit(working_wave,x0, y0, y_norm, x_idx, line_type)
     
        #Calculating Flux: Signal Line Fit

        if type(o1) != type(None):
            dx = x0[2]-x0[1]  #0.16866575
            if line_type == 'Single':
                x_sigsnip = np.where((np.abs((x0-working_wave))/o1[1])<=2.5)[0]
                gauss0=gauss(x0,*o1)
            
            if line_type == 'Balmer':
                x_sigsnip = np.where(np.abs((x0-working_wave))/o1[1]<=2.5 )[0] 
                gauss0 = double_gauss(x0, *o1)
                o1_neg = [o1[0], o1[4], o1[5], o1[3]]
                neg0   = gauss(x0, *o1_neg)
                gauss0_diff = gauss0 - neg0
                y_norm_diff = y_norm[x_sigsnip]-neg0[x_sigsnip]

            if line_type == 'Oxy2':
                con1 = 3728.91/3726.16
                x_sigsnip = np.where(((x0-working_wave)/o1[1]>=-2.5) & ((x0-working_wave*con1)/o1[4]<=2.5))[0]
                #print len(x_sigsnip) #, x_sigsnip
                gauss0 =  oxy2_gauss(x0, *o1)
            
            if line_type == 'Single' or line_type == 'Oxy2':
                flux_g = np.sum((gauss0-o1[3])*dx)    #flux from gaussian distribution 
                flux_s = np.sum((y_norm[x_sigsnip]-o1[3])*dx)  #flux from snipping method (spectral flux)where snip off sigma >2.5

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff*dx)
                flux_s = np.sum(y_norm_diff*dx)


            #Calculating RMS
            ini_sig1, RMS_pix= rms_func(wave,dispersion,lineflag,working_wave, y0, o1[1], scalefact)

            #Line Flag Checking Plots
            pdf= line_flag_check(dataset, fitspath,working_wave, lineflag, wave, y_norm, Spect_1D, line_name,row,col,fig,ax_arr)
            #if rr == 0: print 'o1', o1, flux_g, flux_s, x_sigsnip

            
            #Filling In Arrays
            # print flux_g, type(flux_g)
            flux_g_array[rr] = flux_g 
            flux_s_array[rr] = flux_s
            sigma_array[rr]= o1[1]
            median_array[rr] = o1[3]
            norm_array[rr] = max0
            RMS_array[rr] = ini_sig1
            SN_array[rr] = (flux_s/ini_sig1)
            #if dataset == 'Grid': 
                
                #N_gal_array[rr] = narr_gr[rr]    # asc_tab['area'][rr]
                #R_23_array[rr] = rarr_gr[rr]            #asc_tab['xBar'][rr]
                #O_32_array[rr] = oarr_gr[rr]     #asc_tab['yBar'][rr]
            #if dataset != 'Grid':
            N_gal_array[rr] = asc_tab['area'][rr]
            R_23_array[rr] = asc_tab['xBar'][rr]
            O_32_array[rr] = asc_tab['yBar'][rr]
                
            #R_23_array[rr] = asc_tab['xBar'][rr]
            #O_32_array[rr] = asc_tab['yBar'][rr]

            #print RMS_array[rr]
            #Residuals
            resid = y_norm[x_sigsnip]-gauss0[x_sigsnip] + o1[3]  
            #print len(resid), len(x_sigsnip), len(gauss0)

            

            '''#Error Propogation 
            flux_g_err[rr] = error_prop_chuncodes(flux_g,1)
            flux_s_err[rr] = error_prop_chuncodes(flux_s,1)
            sigma_err[rr] =  error_prop_chuncodes(o1[1],1)
            median_err[rr] = error_prop_chuncodes(o1[3],1)
            norm_err[rr] =   error_prop_chuncodes(max0,1)
            RMS_err[rr] =    error_prop_chuncodes(ini_sig1,1)
            SN_err[rr] =     error_prop_chuncodes(flux_s/ini_sig1,1)'''

        #Plotting
        #if not exists(fitspath+ 'Stacking_Voronoi_Zoomed_Gauss_generalexperiment.pdf'):
       
            emis= t_ax.plot(wave, y_norm,'k', linewidth=0.6, label= 'Emission')
            #if keyword == 'Oxy2':
            t_ax.set_xlim([x1+45,x2-45])
            '''else:
            t_ax.set_xlim([x1,x2])'''


            t_ax.plot(x0,gauss0, 'b--', linewidth= 0.5, label= 'Gauss Fit')
            t_ax.plot(x0[x_sigsnip],resid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')
            t_ax.legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)

            #if keyword == 'Single' or keyword == 'Balmer' : t_ax.axhline(y=o1[3], color='k--', linewidth=0.01)
 
        
            if dataset == 'Grid':
                #txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
                txt0  = 'R_23: %.3f O_32: %.3f\n' % (R_23_array[rr],O_32_array[rr])  #$\overline{x}$:$\overline{y}$:
                txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %.3f\n' % (ini_sig1, RMS_pix, N_gal_array[rr]) 
                txt0 += 'Median: %.3f Sigma: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                txt0 += 'Flux_G: %.3f Flux_S: %.3f' %(flux_g, flux_s) + '\n'
                txt0 += 'S/N: %.3f' %(SN_array[rr])

            else:
                txt0 = r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
                txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])  #$\overline{x}$:$\overline{y}$:
                txt0 += 'RMS: %.3f RMS/pix: %.3f, Scale: %.3f, N: %.3f\n' % (ini_sig1, RMS_pix, asc_tab['scale'][rr], asc_tab['area'][rr]) 
                txt0 += 'Median: %.3f Sigma: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                txt0 += 'Flux_G: %.3f Flux_S: %.3f' %(flux_g, flux_s) + '\n'
                txt0 += 'S/N: %.3f' %(SN_array[rr])
        
       
            t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
            for x in  lambda0: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

            if col == 0:
                t_ax.set_ylabel('Spect_1D')
                t_ax.set_yticklabels([])  #sets y-tick labels 

            if rr == Spect_1D.shape[0]-1 and rr % (nrows*ncols) != nrows*ncols-1:
                for jj in range(col+1, ncols): ax_arr[row,jj].axis('off')
            for kk in range(row+1, nrows):
                for zz in range(ncols): ax_arr[kk,zz].axis('off')
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
                
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == Spect_1D.shape[0]-1: 
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97, hspace=0.05)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(8,8)
            #plt.draw()
            fig.savefig(pdf_pages, format='pdf')

       

        
        #plt.draw()
        #fig.savefig(pdfpages2, format='pdf')
    #endfor
    '''
    #Error Propogation 
    flux_g_err = error_prop_chuncodes(flux_g_array,1)
    flux_s_err = error_prop_chuncodes(flux_s_array,1)
    sigma_err =  error_prop_chuncodes(sigma_array,1)
    median_err = error_prop_chuncodes(median_array,1)
    norm_err =   error_prop_chuncodes(norm_array,1)
    RMS_err =    error_prop_chuncodes(RMS_array,1)
    SN_err =     error_prop_chuncodes(SN_array,1)'''
     
    #Writing Ascii Tables and Fits Tables
    out_ascii = fitspath+'/'+dataset+'_flux_gaussian_'+str(np.int(working_wave))+'.tbl'
    out_fits = fitspath+'/'+dataset+'_Flux_gaussian_'+line_name+'.fits'
    #if not exists(out_ascii):
    n=  ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS', 'S/N')
    n = tuple([line_name + '_' + val for val in n])
    tab0 = Table([ flux_g_array, flux_s_array,sigma_array,median_array,norm_array,RMS_array, SN_array], names=n)
    asc.write(tab0, out_ascii, format='fixed_width_two_line')
        #fits.writeto(out_fits, tab0)
         
    out_ascii_single = fitspath+'/'+dataset+'_Average_R23_O32_Values.tbl'
    if not exists(out_ascii_single):
        n2= ('R_23_Average', 'O_32_Average', 'N_Galaxies', 'RMS')
        tab1 = Table([R_23_array, O_32_array, N_gal_array, RMS_array], names=n2)
        asc.write(tab1, out_ascii_single, format='fixed_width_two_line')



       
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

        
    pdf_pages.close()

    
    print 'Done!'



def zm_general(dataset, fitspath, stack2D, header, wave, lineflag, Spect_1D, dispersion, lambda0, line_type, line_name, s,a,c,s1,a1,s2,a2,tab):
    

    for ii in range(len(lambda0)):
        #print lambda0[ii], line_type[ii], line_name[ii]
        if line_type[ii] == 'Single':
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, Spect_1D,lambda0[ii], lambda0, lineflag, line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

            #print outpdf 'created'
            
        if line_type[ii] == 'Balmer': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, Spect_1D, lambda0[ii], lambda0, lineflag, line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

            #print outpdf 'created'
            
        if line_type[ii] == 'Oxy2': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, Spect_1D, lambda0[ii], lambda0, lineflag, line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

            #print outpdf 'created'
            

        

def line_flag_check(dataset, fitspath,working_wave, lineflag, wave, y_norm, Spect_1D, line_name,row,col,fig,ax_arr):
    #print 'run line_flag_check'

    #New plots for lineflagging
    pdfpages2 = PdfPages(fitspath + '/' +dataset+'_lineflag_check_'+line_name+'.pdf')
    t_ax2 = ax_arr[row,col]
    t_ax2.plot(wave, y_norm, 'k', linewidth=0.6, label='Emission')
    t_ax2.set_xlim([working_wave+150,working_wave-45])
    t_ax2.plot(wave,lineflag,'g', linestyle='dashed', label = 'Lineflag')
    #txt2 =  r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr])
    #t_ax2.annotation(txt2)
    fig.set_size_inches(8,8)
    fig.savefig(pdfpages2, format='pdf')
    pdfpages2.close()



#erg/s/sm^2/A
#plot residuals 



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
