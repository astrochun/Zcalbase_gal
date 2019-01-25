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
RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'

#If 'x0' is infeasible error occurs, check the para_bound values to make sure the expected values are within the range set up upper and lower limits. 

###NOT IN USE###
'''#eventually I need to incoorparate these values into the correct function but right now they are global variables 
s=1.0
a= 1.0
c = 2.0
s1= 1.3
a1= 2.1
s2 = 1.5   ###Changed this value in order to better fit the balmer curves### used to be s2 =1
a2 = -1.0

Works for Single Grids: s=1.0
a= 1.0
c = 2.0
s1= 1.3
a1= 1.5
s2 = 1.5   ###Changed this value in order to better fit the balmer curves### used to be s2 =1
a2 = 1.8'''
###   

def line_flag_check(dataset, fitspath,working_wave, lineflag, wave, y_norm, Spect_1D, line_name,row,col,fig,ax_arr):
    #print 'Under Construction '

    #New plots for lineflagging
    pdfpages2 = PdfPages(fitspath + '/' +dataset+'_lineflag_check_'+line_name+'.pdf')
    t_ax2 = ax_arr[row,col]
    t_ax2.plot(wave, y_norm, 'k', linewidth=0.4, label='Emission')
    t_ax2.set_xlim([working_wave+150,working_wave-45])
    t_ax2.plot(wave,lineflag,'g', linestyle='dashed', label = 'Lineflag')
    #txt2 =  r'xnode=%.3f  ynode=%.3f' % (asc_tab['xnode'][rr], asc_tab['ynode'][rr])
    #t_ax2.annotation(txt2)
    fig.set_size_inches(8,8)
    fig.savefig(pdfpages2, format='pdf')
    pdfpages2.close()
    

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

def get_gaussian_fit(dataset, s2, working_wave,x0, y0, y_norm, x_idx, RMS, line_type):

    med0 = np.median(y_norm[x_idx])
    max0 = np.max(y_norm[x_idx]) - med0
    sigma = np.repeat(RMS,len(x0))
    #print sigma
    
    fail = 0


    ###Single Emission Line###
    if line_type == 'Single': 
        p0 = [working_wave, 1.0, max0, med0] #must have some reasonable values
       
        para_bound = ((working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0)),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0)))

        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma= sigma[x_idx], bounds = para_bound)
        except ValueError:
            print 'fail'
            fail = 1

    
    ###Double Balmer Emission Line###
    ###initial para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0,0)###
    if line_type == 'Balmer': 
        '''p0 = [working_wave, 1.0, max0, med0, s2, -0.05*max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0),0.0, -max0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),25.0,0.0)'''

        if dataset == 'R23_Grid' or dataset =='O32_Grid':    
            p0 = [working_wave, 1.0, max0, med0, s2, -0.05*max0] #must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0),0.0, -max0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),30.0,0.0)

        if dataset == 'Grid':
            #print 'med0:', med0, 'max0:',max0
            p0 = [working_wave, 1.0, max0, med0, s2, -0.25*med0] #must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0),0.0, -max0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),30.0,0.0)

        if dataset == 'Voronoi10' or dataset =='Voronoi14' or dataset== 'Voronoi20':
            #print 'med0:', med0, 'max0:', max0
            p0 = [working_wave, 1.0, max0, med0, s2, -0.50*max0] #must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0),0.0, -max0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),30.0,0) 

        #print para_bound
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=sigma[x_idx], bounds = para_bound)
        except ValueError:
            print 'fail'
            fail = 1


    ###OxygenII Emission Line###
    if line_type == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75*max0, med0, 1.0, max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, 0.0),(working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0, 100.0)

        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx], p0=p0, sigma=sigma[x_idx], bounds = para_bound)
        except ValueError:
            print 'fail'
            fail = 1
   
    if not fail:
        #print 'o1:', o1
        return o1, med0, max0
    else:
        return None, med0, max0


def rms_func(wave,dispersion,lineflag, lambda_in, y0, scalefact, sigma_array):
    x_idx = np.where((np.abs(wave-lambda_in)<=100) & (lineflag==0))[0]
    sigma = np.std(y0[x_idx])
    RMS_pix = sigma*dispersion/scalefact
    if sigma_array == 0:
        RMS = sigma/scalefact
        return RMS
    else: 
        pix =  5* sigma_array / dispersion 
        s_pix = np.sqrt(pix)
        ini_sig= s_pix * sigma * dispersion
        return ini_sig/scalefact, RMS_pix


def error_prop_chuncodes(x_arr, dx):
    #x_arr = np.random.random_integers(10,10)
    x_pdf = random_pdf(x_arr,0.5)
    return x_pdf
    
#for each individual stack
#electron temperature and the R23 and O32 values 
#Plotting Zoomed 
def zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, working_wave, lambda0, lineflag,  y_correction='', line_type = '',outpdf='', line_name=''):
    
    asc_tab = asc.read(tab)
    pdf_pages = PdfPages(outpdf)
    nrows = 4
    ncols = 4
    x_idx = np.where((wave>=(working_wave-100)) & (wave<=(working_wave+100)))[0] 
    x0 = wave#[x_idx]
    scalefact = 1e-17

    #Initializing Arrays

    flux_g_array = np.zeros(stack2D.shape[0])
    flux_s_array = np.zeros(stack2D.shape[0])
    sigma_array = np.zeros(stack2D.shape[0])
    median_array = np.zeros(stack2D.shape[0])
    norm_array = np.zeros(stack2D.shape[0])
    RMS_array = np.zeros(stack2D.shape[0])
    SN_array = np.zeros(stack2D.shape[0])
    N_gal_array = np.zeros(stack2D.shape[0])
    R_23_array = np.zeros(stack2D.shape[0])
    O_32_array = np.zeros(stack2D.shape[0])

    for rr in range(stack2D.shape[0]):
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

        y_smooth = movingaverage_box1D(y_norm, 4, boundary='extend')
      
        RMS= rms_func(wave,dispersion,lineflag,working_wave, y0, scalefact, 0)
        
        if y_correction == '': o1, med0, max0  = get_gaussian_fit(dataset, s2, working_wave,x0, y0, y_norm, x_idx, RMS, line_type)

        if y_correction == 'y_smooth': o1, med0, max0 = get_gaussian_fit(dataset, s2, working_wave,x0, y0, y_smooth, x_idx, RMS, line_type)
        
        #Calculating Flux: Signal Line Fit

        if type(o1) != type(None):
            dx = x0[2]-x0[1] 
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
                gauss0 =  oxy2_gauss(x0, *o1)
            
            if line_type == 'Single' or line_type == 'Oxy2':
                flux_g = np.sum((gauss0-o1[3])*dx)    #flux from gaussian distribution 
                flux_s = np.sum((y_norm[x_sigsnip]-o1[3])*dx)  #flux from snipping method (spectral flux)where snip off sigma >2.5

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff*dx)
                flux_s = np.sum(y_norm_diff*dx)


            #Calculating RMS
            ini_sig1, RMS_pix= rms_func(wave,dispersion,lineflag,working_wave, y0, scalefact, o1[1])

            #Line Flag Checking Plots
            #line_flag_check(dataset, fitspath,working_wave, lineflag, wave, y_norm, stack2D, line_name,row,col,fig,ax_arr)
           

            
            #Filling In Arrays
            # print flux_g, type(flux_g)
            flux_g_array[rr] = flux_g 
            flux_s_array[rr] = flux_s
            sigma_array[rr]= o1[1]
            median_array[rr] = o1[3]
            norm_array[rr] = max0
            RMS_array[rr] = ini_sig1
            SN_array[rr] = (flux_s/ini_sig1)
          
            #if dataset != 'Grid':
            N_gal_array[rr] = asc_tab['area'][rr]
            R_23_array[rr] = asc_tab['xBar'][rr]
            O_32_array[rr] = asc_tab['yBar'][rr]
        
            #Residuals
            resid = y_norm[x_sigsnip]-gauss0[x_sigsnip] + o1[3]  

            

        #Plotting
       
            emis= t_ax.plot(wave, y_smooth,'k', linewidth=0.4, label= 'Emission')
            if y_correction =='y_smooth': t_ax.plot(x0,gauss0, 'm--', linewidth= 0.5, label= 'Gauss Fit')
            else: t_ax.plot(x0,gauss0, 'b--', linewidth= 0.5, label= 'Gauss Fit')
            t_ax.plot(x0[x_sigsnip],resid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')
            t_ax.plot(wave,lineflag,'g', linestyle='dashed', linewidth = 0.1, label = 'Lineflag')
            #t_ax.legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2, fontsize = 3)
            t_ax.set_xlim(x1,x2)
 
            
            
            
            if dataset == 'Grid' or dataset=='O32_Grid' or dataset =='R23_Grid':
                if line_type == 'Balmer': 
                    txt0  = 'ID: %i, R_23: %.3f O_32: %.3f\n' % (asc_tab['ID'][rr], R_23_array[rr],O_32_array[rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %.3f\n' % (ini_sig1, RMS_pix, N_gal_array[rr]) 
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

                if line_type == 'Single' or line_type =='Oxy2': 
                    txt0  = 'ID: %i, R_23: %.3f O_32: %.3f\n' % (asc_tab['ID'][rr], R_23_array[rr],O_32_array[rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %i\n' % (ini_sig1, RMS_pix, N_gal_array[rr]) 
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f'% (o1[3], o1[1], max0, o1[2]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

            else:
                if line_type =='Balmer':
                    txt0 = r'ID: %i  xnode=%.3f  ynode=%.3f' % (asc_tab['ID'][rr], asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
                    txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, Scale: %.3f, N: %.3f\n' % (ini_sig1, RMS_pix, asc_tab['scale'][rr], N_gal_array[rr]) 
                    txt0 += 'Median: %.3f $\sigma$: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

                if line_type =='Single' or line_type=='Oxy2':
                    txt0 = r'ID: %i  xnode=%.3f  ynode=%.3f' % (asc_tab['ID'][rr], asc_tab['xnode'][rr], asc_tab['ynode'][rr]) + '\n'
                    txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['xBar'][rr], asc_tab['yBar'][rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, Scale: %.3f, N: %i \n' % (ini_sig1, RMS_pix, asc_tab['scale'][rr], N_gal_array[rr]) 
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f'% (o1[3], o1[1], max0,o1[2]) + '\n'
                    txt0 += 'Flux_G: %.3f Flux_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])
        
       
            t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
            for x in  lambda0: t_ax.axvline(x=x, linewidth= 0.3, color= 'k')

            if col == 0:
                t_ax.set_ylabel('Spect_1D')
            else: t_ax.set_yticklabels([])  #sets y-tick labels 

            if row != nrows-1 and rr != stack2D.shape[0]-1:
                t_ax.set_xticklabels([])
                
            #if rr == Spect_1D.shape[0]-1 and rr % (nrows*ncols) != nrows*ncols-1:
            #    for jj in range(col+1, ncols): ax_arr[row,jj].axis('off')
            #for kk in range(row+1, nrows):
            #    for zz in range(ncols): ax_arr[kk,zz].axis('off')
            
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == stack2D.shape[0]-1: 
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97, hspace=0.05)
            
            #if t_ax > Spect_1D.shape[0]: fig.delaxes(t_ax) #use a for loop 

            fig.set_size_inches(8,8)
            #plt.draw()
            fig.savefig(pdf_pages, format='pdf')

       

        
        #plt.draw()
        #fig.savefig(pdfpages2, format='pdf')
    #endfor
    
     
    #Writing Ascii Tables and Fits Tables
    ID = asc_tab['ID'].data
    out_ascii = fitspath+'/'+dataset+'_flux_gaussian_'+str(np.int(working_wave))+'.tbl'
    out_fits = fitspath+'/'+dataset+'_Flux_gaussian_'+line_name+'.fits'
    #if not exists(out_ascii):
    n=  ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS', 'S/N')
    n = tuple([line_name + '_' + val for val in n])
    tab0 = Table([flux_g_array, flux_s_array,sigma_array,median_array,norm_array,RMS_array, SN_array], names=n)
    asc.write(tab0, out_ascii, format='fixed_width_two_line')
        #fits.writeto(out_fits, tab0)
         
    out_ascii_single = fitspath+'/'+dataset+'_Average_R23_O32_Values.tbl'
    #if not exists(out_ascii_single):
    n2= ('ID','R_23_Average', 'O_32_Average', 'N_Galaxies', 'RMS')
    tab1 = Table([ID, R_23_array, O_32_array, N_gal_array, RMS_array], names=n2)
    asc.write(tab1, out_ascii_single, format='fixed_width_two_line')

        
    pdf_pages.close()

    
    print 'Done!'

    fig.clear()

def zm_general(dataset, fitspath, stack2D, header, wave, lineflag, dispersion, lambda0, line_type, line_name,  y_correction,s,a,c,s1,a1,s2,a2,tab):
    

    for ii in range(len(lambda0)):
        #Single Gaussian Fit 
        if line_type[ii] == 'Single':
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave,lambda0[ii], lambda0, lineflag, y_correction='', line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

        #Balmer Line Fit  
        if line_type[ii] == 'Balmer': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, lambda0[ii], lambda0, lineflag, y_correction=y_correction, line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

         
        #Oxy2 Line Fit     
        if line_type[ii] == 'Oxy2': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print outpdf
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, header, dispersion,  s,a,c,s1,a1,s2,a2, wave, lambda0[ii], lambda0, lineflag, y_correction='', line_type= line_type[ii], outpdf=outpdf, line_name=line_name[ii])

            

    




