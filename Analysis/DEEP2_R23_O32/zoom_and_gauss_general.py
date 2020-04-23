
"""
PERFORMS THE GAUSSIAN FITTING ON ALL EMISSIONS LINES
CALLED EVERY TIME THE GENERAL FUNCTIONS ARE RUN
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column
from pylab import subplots_adjust
from scipy.optimize import curve_fit

# Import error propagation codes from chun_codes
from chun_codes import random_pdf

from Metallicity_Stack_Commons.Metallicity_Stack_Commons.analysis.fitting import gauss, double_gauss, oxy2_gauss
from Metallicity_Stack_Commons.Metallicity_Stack_Commons.analysis.fitting import movingaverage_box1D, rms_func
from Metallicity_Stack_Commons.Metallicity_Stack_Commons import lambda0, line_name, line_type

'''
Debugging Note:
   If 'x0' is infeasible error occurs, check the para_bound values to
   make sure the expected values are within the range set up upper and
   lower limits.
'''


def line_flag_check(dataset, fitspath, working_wave, lineflag, wave, y_norm,
                    line_name0, row, col, fig, ax_arr):

    # New plots for line flagging

    out_pdf = '%s/%s_lineflag_check_%s.pdf' % (fitspath, dataset, line_name0)
    pdf_pages2 = PdfPages(out_pdf)

    t_ax2 = ax_arr[row, col]
    t_ax2.plot(wave, y_norm, 'k', linewidth=0.4, label='Emission')
    t_ax2.set_xlim([working_wave+150, working_wave-45])
    t_ax2.plot(wave, lineflag, 'g', linestyle='dashed', label='Lineflag')

    fig.set_size_inches(8, 8)
    fig.savefig(pdf_pages2, format='pdf')
    pdf_pages2.close()


def get_gaussian_fit(dataset, s2, working_wave, x0, y_norm, x_idx, RMS, line_type0):

    med0 = np.median(y_norm[x_idx])
    max0 = np.max(y_norm[x_idx]) - med0
    sigma = np.repeat(RMS, len(x0))

    fail = 0

    # curve_fit(gauss function, x0 ,y0, initial_guesses, sigma, parameter_bounds)
    # Single Emission Line
    if line_type0 == 'Single':
        p0 = [working_wave, 1.0, max0, med0]  # must have some reasonable values

        para_bound = ((working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0)),
                      (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0)))

        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0,
                               sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1

    # Double Balmer Emission Line
    '''
    initial para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),
                         (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),10.0,0)
    '''
    if line_type0 == 'Balmer':
        '''
        p0 = [working_wave, 1.0, max0, med0, s2, -0.05*max0] #must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0),0.0, -max0),
                     (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0),25.0,0.0)
        '''

        if dataset == 'R23_Grid' or dataset == 'O32_Grid':
            p0 = [working_wave, 1.0, max0, med0, s2, -0.05*max0]  # must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),\
                         (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0), 30.0, 0.0)

        if dataset == 'Grid':
            p0 = [working_wave, 1.0, max0, med0, s2, -0.25*med0]  # must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),\
                         (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0), 30.0, 0.0)

        if dataset == 'Voronoi10' or dataset == 'Voronoi14' or dataset == 'Voronoi20' \
                or dataset == 'Double_Bin' or dataset == 'n_Bins':

            p0 = [working_wave, 1.0, max0, med0, s2, -0.5*med0]  # must have some reasonable values
            para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, -med0),\
                         (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0), 30.0, 0)

        # Attempt fit
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1

    # OxygenII Emission Line
    if line_type0 == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75*max0, med0, 1.0, max0]  # must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, 0.0),\
                     (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0), 10.0, 100.0)

        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            print('fail')
            fail = 1
   
    if not fail:
        return o1, med0, max0
    else:
        return None, med0, max0


def equi_width_func(pos_comp, neg0, gauss0, x0, wave, y_norm):
    """
    Purpose:
      Equivalent width correction/computation


    Notes:
      bottom side of the iceburg / continuum
      take negative component of gauss and subtract off the positive component

      total gauss - the positive component?
      x = double gaussian of o1[0,1,2] = 0
      x-o1[3] from [-2.5sigma to 2.5 sigma]
      equivalent width in terms of Angstroms
      update plots
    """

    plt.plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
    plt.plot(x0, neg0, 'b', linewidth=0.25, label='Negative Component')
    plt.plot(x0, pos_comp, 'r', linewidth=0.25, label='Positive Component')
    plt.plot(x0, gauss0, 'g', linewidth=0.25, label='Gauss Fit')


# For each individual stack
# Electron temperature and the R23 and O32 values
# Plotting Zoomed
def zoom_gauss_plot(dataset, fitspath, tab, stack2D, dispersion, s2, wave,
                    working_wave, lineflag,  y_correction='', line_type='',
                    outpdf='', line_name=''):

    asc_tab = asc.read(tab)
    pdf_pages = PdfPages(outpdf)
    nrows = 4
    ncols = 4
    x_idx = np.where((wave >= (working_wave-100)) & (wave <= (working_wave+100)))[0]
    x0 = wave
    scalefact = 1e-17

    ID = asc_tab['bin_ID'].data

    # Initializing Arrays
    flux_g_array = np.zeros(stack2D.shape[0])
    flux_s_array = np.zeros(stack2D.shape[0])
    flux_neg_array = np.zeros(stack2D.shape[0])
    sigma_array = np.zeros(stack2D.shape[0])
    median_array = np.zeros(stack2D.shape[0])
    norm_array = np.zeros(stack2D.shape[0])
    RMS_array = np.zeros(stack2D.shape[0])
    SN_array = np.zeros(stack2D.shape[0])
    N_gal_array = np.zeros(stack2D.shape[0])
    R_23_array = np.zeros(stack2D.shape[0])
    O_32_array = np.zeros(stack2D.shape[0])

    #Initializing Arrays for Balmer Graphing
    xbar_array = np.zeros(stack2D.shape[0])
    sig1_array = np.zeros(stack2D.shape[0])
    pos_amp_array = np.zeros(stack2D.shape[0])
    const_array = np.zeros(stack2D.shape[0])
    sig2_array = np.zeros(stack2D.shape[0])
    neg_amp_array = np.zeros(stack2D.shape[0])

    for rr in range(stack2D.shape[0]):
        y0 = stack2D[rr]
        y_norm = y0/scalefact

        row = rr // nrows % ncols
        col = rr % ncols
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False)

        t_ax = ax_arr[row, col]
        
        x1 = working_wave-100
        x2  = working_wave+100

        y_smooth = movingaverage_box1D(y_norm, 4, boundary='extend')

        RMS_ang = rms_func(wave, dispersion, working_wave, y0, 0, lineflag)

        if y_correction == '':
            o1, med0, max0  = get_gaussian_fit(dataset, s2, working_wave, x0,
                                               y_norm, x_idx, RMS_ang,
                                               line_type)

        if y_correction == 'y_smooth':
            o1, med0, max0 = get_gaussian_fit(dataset, s2, working_wave,x0,
                                              y_smooth, x_idx, RMS_ang,
                                              line_type)

        # Calculating Flux: Signal Line Fit
        if type(o1) != type(None):
            dx = x0[2]-x0[1] 
            if line_type == 'Single':
                x_sigsnip = np.where((np.abs((x0-working_wave))/o1[1])<=2.5)[0]
                gauss0=gauss(x0, *o1)
            
            if line_type == 'Balmer':
                x_sigsnip = np.where(np.abs((x0-working_wave))/o1[1]<=2.5 )[0]
                gauss0 = double_gauss(x0, *o1)

                o1_neg = [o1[0], o1[4], o1[5], o1[3]]
                neg0   = gauss(x0, *o1_neg)
                gauss0_diff = gauss0 - neg0
                y_norm_diff = y_norm[x_sigsnip]-neg0[x_sigsnip]

                # Equivalent Width calculation
                pos_comp = gauss0-neg0
                neg_comp = neg0[x_sigsnip]
                flux_neg = np.sum(neg_comp*dx)
                equiv_wid = flux_neg/o1[3]

                # Keeping this here for possible later use. Might belong in MSC balmer module
                # pdfpages3 = PdfPages(fitspath + '/' +dataset+'_PlottingBalmer_'+line_name+'.pdf')
                # equi_width_func(fitspath, dataset, working_wave, pos_comp, neg0, gauss0, x0, wave, y_norm, line_type,pdfpages3)

            if line_type == 'Oxy2':
                con1 = 3728.91/3726.16
                x_sigsnip = np.where(((x0-working_wave)/o1[1]>=-2.5) &
                                     ((x0-working_wave*con1)/o1[4]<=2.5))[0]
                gauss0 =  oxy2_gauss(x0, *o1)
            
            if line_type == 'Single' or line_type == 'Oxy2':
                flux_g = np.sum((gauss0-o1[3])*dx)            # Flux from gaussian distribution
                flux_s = np.sum((y_norm[x_sigsnip]-o1[3])*dx) # Flux from snipping method (spectral flux)where snip off sigma >2.5

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff*dx)
                flux_s = np.sum(y_norm_diff*dx)


            #Calculating RMS
            RMS_tot, RMS_pix= rms_func(wave, dispersion,working_wave,y0,o1[1],lineflag)

            #Line Flag Checking Plots
            #line_flag_check(dataset, fitspath,working_wave, lineflag, wave, y_norm, stack2D, line_name,row,col,fig,ax_arr)

            #Filling In Arrays
            # print flux_g, type(flux_g)
            flux_g_array[rr] = flux_g
            flux_s_array[rr] = flux_s
            sigma_array[rr]= o1[1]
            median_array[rr] = o1[3]
            norm_array[rr] = max0
            RMS_array[rr] = RMS_tot
            SN_array[rr] = (flux_s/RMS_tot)
            if line_type == 'Balmer': flux_neg_array[rr] = flux_neg

            N_gal_array[rr] = asc_tab['N_stack'][rr]
            R_23_array[rr] = asc_tab['logR23_avg'][rr]
            O_32_array[rr] = asc_tab['logO32_avg'][rr]
        
            # Residuals
            if line_type == 'Oxy2': 
                x_sigsnip_2 = np.where(np.abs(x0-3727)/o1[1]<=4.0)[0]
                resid = y_norm[x_sigsnip_2]-gauss0[x_sigsnip_2] + o1[3]
            else: 
                x_sigsnip_2 = np.where(np.abs((x0-working_wave))/o1[1]<=3.0 )[0]
                resid = y_norm[x_sigsnip_2]-gauss0[x_sigsnip_2] + o1[3]

            


            # Filling in Balmer Arrays
            if line_type == "Single":
                xbar_array[rr] = o1[0]
                sig1_array[rr] = o1[1]
                pos_amp_array[rr]= o1[2]
                const_array[rr] = o1[3]
                 
            if line_type == 'Balmer' or line_type == 'Oxy2':
                xbar_array[rr] = o1[0]
                sig1_array[rr] = o1[1]
                pos_amp_array[rr]= o1[2]
                const_array[rr] = o1[3]
                sig2_array[rr] = o1[4]
                neg_amp_array[rr] = o1[5]

            # Plotting
            if y_correction == 'y_smooth':
                emis = t_ax.plot(wave, y_smooth,'k', linewidth=0.3, label= 'Emission')
                t_ax.plot(x0, gauss0, 'm', linewidth=0.25, label='Gauss Fit')
            if y_correction == '':
                emis = t_ax.plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
                t_ax.plot(x0,gauss0, 'b', linewidth= 0.25, label= 'Gauss Fit')

            t_ax.plot(x0[x_sigsnip_2],resid, 'r', linestyle='dashed', linewidth=0.2,
                      label= 'Residuals')
            #t_ax.plot(wave,lineflag,'g', linestyle='dashed', linewidth=0.1, label='Lineflag')
            t_ax.set_xlim(x1+50,x2-50)

            if line_name == 'OIII_4363':
                t_ax.set_ylim(0,1)
            
            if dataset == 'Grid' or dataset=='O32_Grid' or dataset =='R23_Grid' \
                    or dataset =='Double_Bin' or dataset =='n_Bins':
                if line_type == 'Balmer': 
                    txt0  = 'Line: %.3f, ID: %i, R_23: %.3f O_32: %.3f\n' % (o1[0], ID[rr], asc_tab['logR23_min'][rr], asc_tab['logO32_min'][rr]) + '\n'
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %.3f\n' % (RMS_tot, RMS_pix, N_gal_array[rr]) 
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

                if line_type == 'Single' or line_type =='Oxy2': 
                    txt0  = 'Line: %.3f, ID: %i, R_23: %.3f O_32: %.3f\n' % (o1[0], ID[rr], asc_tab['logR23_min'][rr], asc_tab['logO32_min'][rr]) + '\n'
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %i\n' % (RMS_tot, RMS_pix, N_gal_array[rr]) 
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f'% (o1[3], o1[1], max0, o1[2]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

            else:
                if line_type =='Balmer':
                    txt0 = r'Line: %.3f, ID: %i  xnode=%.3f  ynode=%.3f' % (o1[0], ID[rr], asc_tab['logR23_min'][rr], asc_tab['logO32_min'][rr]) + '\n'
                    txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['logR23_avg'][rr], asc_tab['logO32_avg'][rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %.3f\n' % (RMS_tot, RMS_pix, N_gal_array[rr])
                    txt0 += 'Median: %.3f $\sigma$: %.3f  Norm: %.3f'% (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) + '\n'
                    txt0 += 'F_G: %.3f F_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

                if line_type =='Single' or line_type=='Oxy2':
                    txt0 = r'Line: %.3f, ID: %i  xnode=%.3f  ynode=%.3f' % (o1[0], ID[rr], asc_tab['logR23_min'][rr], asc_tab['logO32_min'][rr]) + '\n'
                    txt0 += 'R_23: %.3f O_32: %.3f\n' % (asc_tab['logR23_avg'][rr], asc_tab['logO32_avg'][rr])
                    txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %i \n' % (RMS_tot, RMS_pix, N_gal_array[rr])
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f'% (o1[3], o1[1], max0,o1[2]) + '\n'
                    txt0 += 'Flux_G: %.3f Flux_S: %.3f' %(flux_g, flux_s) + '\n'
                    txt0 += 'S/N: %.3f' %(SN_array[rr])

            
            #t_ax.annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right',fontsize= '5')
            for x in  lambda0: t_ax.axvline(x=x, linewidth= 0.1, color= 'k', linestyle = '--')

            txt1 = 'Intensity' #Units: ' + r'$ergs *s^{-1} *cm^{-2} *\AA$' + '\n'
            #txt1 += 'Scale factor = 1e-17'
            
            if col == 0:
                t_ax.set_ylabel(txt1)
            else: t_ax.set_yticklabels([])  # sets y-tick labels

            if row != nrows-1 and rr != stack2D.shape[0]-1:
                t_ax.set_xticklabels([])
        else:
            t_ax.plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
            t_ax.set_xlim(x1+50,x2-50)
            
            
        if (rr % (nrows*ncols) == nrows*ncols-1) or rr == stack2D.shape[0]-1: 
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97, hspace=0.05)

            fig.set_size_inches(8,8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()
    # endfor


    # Writing Ascii Tables and Fits Tables
    # Repetive Columns have been added: The last four ot six columns will be used for individual graphing
    out_ascii = fitspath+'/'+dataset+'_flux_gaussian_'+str(np.int(working_wave))+'.tbl'
    # out_fits = fitspath+'/'+dataset+'_Flux_gaussian_'+line_name+'.fits'


    
    if line_type == 'Single':
        n = ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS',
             'S/N', 'X_bar', 'Pos_Sig', 'Pos_Amp', 'Const')

        n = tuple([line_name + '_' + val for val in n])
        tab0 = Table([flux_g_array, flux_s_array, sigma_array, median_array,
                      norm_array,RMS_array, SN_array, xbar_array, sig1_array,
                      pos_amp_array, const_array], names=n)
        asc.write(tab0, out_ascii, format='fixed_width_two_line')

    if line_type == 'Balmer' or line_type == 'Oxy2': 
         n=  ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS',
              'S/N', 'X_bar', 'Pos_Sig', 'Pos_Amp', 'Const', 'Neg_Sig', 'Neg_Amp')

         n = tuple([line_name + '_' + val for val in n])
         tab0 = Table([flux_g_array, flux_s_array, sigma_array, median_array,
                       norm_array, RMS_array, SN_array, xbar_array, sig1_array,
                       pos_amp_array, const_array, sig2_array, neg_amp_array], names=n)

         if line_type == 'Balmer':
             print('Adding an Equ_Width Column')
             names = 'EW_'+str(np.int(working_wave))+'_abs'
             equ_add = Column(name=names, data=flux_neg_array)
             tab0.add_column(equ_add, 2)
    asc.write(tab0, out_ascii, format='fixed_width_two_line')


    
    out_ascii_single = fitspath+'/'+dataset+'_Average_R23_O32_Values.tbl'


    n2= ('bin_ID','logR23_avg', 'logO32_avg', 'N_stack')
    tab1 = Table([ID, R_23_array, O_32_array, N_gal_array], names=n2)
    asc.write(tab1, out_ascii_single, format='fixed_width_two_line')



                           

    pdf_pages.close()

    # pdfpages3.close()
    
    
    print('Done!')

    fig.clear()
def zm_general(dataset, fitspath, stack2D, wave, lineflag, dispersion, y_correction,s,a,c,s1,a1,s2,a2,tab):


    for ii in range(len(lambda0)):
        # Single Gaussian Fit
        if line_type[ii] == 'Single':
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print(outpdf)
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, dispersion, s2,
                            wave, lambda0[ii], lineflag, 
                            y_correction=y_correction, line_type=line_type[ii],
                            outpdf=outpdf, line_name=line_name[ii])

        # Balmer Line Fit
        if line_type[ii] == 'Balmer': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print(outpdf)
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, dispersion, s2,
                            wave, lambda0[ii], lineflag,
                            y_correction=y_correction, line_type=line_type[ii],
                            outpdf=outpdf, line_name=line_name[ii])

         
        # Oxy2 Line Fit
        if line_type[ii] == 'Oxy2': 
            outpdf = fitspath+dataset+'_Zoomed_Gauss_'+line_name[ii]+'.pdf'
            print(outpdf)
            zoom_gauss_plot(dataset, fitspath, tab, stack2D, dispersion, s2,
                            wave, lambda0[ii], lineflag,
                            y_correction=y_correction, line_type=line_type[ii],
                            outpdf=outpdf, line_name=line_name[ii])
