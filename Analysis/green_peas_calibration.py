"""
green_peas_calibration
====

Plots R23 and O32 line ratios and compare it to Jiang et al. (2018)
calibration that is based on green pea galaxies
"""

import sys, os

from chun_codes import systime, match_nosort_str

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

import matplotlib.pyplot as plt
import glob

from scipy.optimize import curve_fit

from scipy.interpolate import interp1d

from astropy.table import Table
from astropy import log

from astropy.io import fits

jiang18_coeffs = [-24.135, 6.1523, -0.37866, -0.147, -7.071]

def O32_OH_fit(xy, a, b, c, d, e):
    '''
    Main functional code that determine log(R23) from log(O32) and 12+log(O/H)

    Parameters
    ----------
     x : 12+log(O/H)
     y : log([OIII]/[OII])
    '''

    x = xy[0]
    y = xy[1]

    logR23 = a + b*x + c * x**2 - d * (e + x) * y

    return logR23

def jiang18(x, y):
    '''
    Function to return log(R23) based on metallicity and [OIII]/[OII] flux ratio

    Parameters
    ----------
     x : 12+log(O/H)
     y : log([OIII]/[OII])
    '''

    logR23 = O32_OH_fit(xy, *jiang18_coeffs)

    return logR23
#enddef

def plot_differences(lR23, lO32, OH, lO32_all, out_diff_pdf, bin_start, bin_end, n_bins=4,
                     lR23_err=[], OH_err=[], OH_range=[], dR23_range=[], marker=[], label=[]):
    '''
    Plot differences between Jiang18 R23 vs observed R23 as a function of metallicity
    '''

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    #Plotting
    ctype = ['red','magenta','green','cyan','blue','black']

    if len(marker) == 0:
        marker = ['o'] * n_sample

    diff0 = []
    for nn in range(n_sample):
        jiang_R23 = O32_OH_fit((OH[nn], lO32[nn]), *jiang18_coeffs)

        # Label in upper left the points
        if len(label) != 0:
            x1 = OH_range[0] + 0.025*(OH_range[1]-OH_range[0])
            y1 = dR23_range[1] - (nn*0.035 + 0.05)*(dR23_range[1]-dR23_range[0])
            x2 = OH_range[0] + 0.035*(OH_range[1]-OH_range[0])
            y2 = dR23_range[1] - (nn*0.035 + 0.0525)*(dR23_range[1]-dR23_range[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1],[y1], marker=marker[nn], color='black')

        for ii in range(n_bins):
            y_ii_min = bin_start[ii]
            y_ii_max = bin_end[ii]
            idx = np.where((lO32[nn] >= y_ii_min) & (lO32[nn] <= y_ii_max))[0]

            ii_label = ''
            if nn == 0: #n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) & (lO32_all <= y_ii_max))[0]
                ii_label = r' %.2f < $\log(O_{32})$ < %.2f, N = %i' % (y_ii_min, y_ii_max,
                                                                       len(idx_all))

            if len(idx) > 0:
                i_diff = lR23[nn][idx] - jiang_R23[idx]
                ax.scatter(OH[nn][idx], i_diff, color=ctype[ii], marker=marker[nn],
                           edgecolor='none', alpha=0.5, label=ii_label)

                diff0 += list(lR23[nn][idx] - jiang_R23[idx])

                if len(OH_err) != 0:
                    ax.errorbar(OH[nn][idx], i_diff, xerr=OH_err[nn][:,idx], mfc='none',
                                ecolor=ctype[ii], capsize=0, alpha=0.25, fmt=None, label=None)

                if len(lR23_err) != 0:
                    ax.errorbar(OH[nn][idx], i_diff, yerr=lR23_err[nn][:,idx], mfc='none',
                                ecolor=ctype[ii], capsize=0, alpha=0.25, fmt=None, label=None)

    # Draw horizontal line at zero:
    ax.axhline(y=0, c='k', linestyle='dashed')

    # Compute statistics
    med0 = np.median(diff0)
    avg0 = np.average(diff0)
    sig0 = np.std(diff0)
    ax.axhline(y=avg0, c='r', linestyle='dotted')
    ax.axhline(y=med0, c='b', linestyle='dotted')

    an_txt  = r'$<\Delta_{R_{23}}>$ : %0.2f' % avg0 + '\n'
    an_txt += r'$\tilde\Delta_{R_{23}}$ : %0.2f' % med0 + '\n'
    an_txt += r'$\sigma$ : %0.2f' % sig0
    ax.annotate(an_txt, [0.155,0.015], xycoords='axes fraction', va='bottom', ha='right',
                fontsize=10)

    if len(OH_range)   != 0: ax.set_xlim(OH_range)
    if len(dR23_range) != 0: ax.set_ylim(dR23_range)

    ax.set_xlabel(r'$12+\log({\rm O/H})_{T_e}$')
    ax.set_ylabel(r'$\Delta_{R_{23}} \equiv \log(R_{23}) - \log(R_{23})_{\rm J18}$')
    ax.minorticks_on()
    leg = ax.legend(loc='upper right', scatterpoints=1, fontsize=8, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.08, top=0.97)

    fig.savefig(out_diff_pdf)
#enddef

def main(lR23, lO32, OH, out_pdf, n_bins=4, lR23_err=[], OH_err=[], xra=[], yra=[],
         marker=[], label=[], dR23_range=[-0.3,0.3], fit=False, silent=False, verbose=True):

    '''
    Main function to plot dataset against Jiang+ (2018) calibration

    Parameters
    ----------
    lR23 : log(R_23)
    lO32 : log([OIII]/[OII])
    OH   : 12+log(O/H)

    out_pdf : full path for output PDF

    fit : boolean
      Turn on fitting. Default: False -> Uses Jiang+2018 relation

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 November 2018
    '''

    if silent == False: log.info('### Begin main : '+systime())

    fig, ax = plt.subplots()

    n_sample = len(lR23)

    min1, max1 = np.zeros(n_sample), np.zeros(n_sample)
    OH_min1, OH_max1 = np.zeros(n_sample), np.zeros(n_sample)
    for nn in range(n_sample):
        min1[nn] = np.min(lO32[nn])
        max1[nn] = np.max(lO32[nn])

        good = np.where(np.isfinite(OH[nn]) == True)[0]
        OH_min1[nn] = np.min(OH[nn][good])
        OH_max1[nn] = np.max(OH[nn][good])

    #print 'min1:', min1,   'max1:', max1
    # Grid of 12+log(O/H)
    if len(yra) == 0:
        x_arr = np.arange(min(OH_min1),max(OH_max1),0.05)
    else:
        x_arr = np.arange(yra[0], yra[1], 0.05)
    #print 'x_arr:', x_arr

    lO32_all = np.array([])
    for nn in range(n_sample):
        lO32_all = np.append(lO32_all, lO32[nn])


    ###Binning based on O32 data because of the way the Jiang calibration works 
    sort0   = np.argsort(lO32_all)
    y_sort0 = lO32_all[sort0]

    bin_pts   = np.int(len(lO32_all)/n_bins)
    r_bin_pts = np.int(np.round(len(lO32_all)/float(n_bins)))

    bin_start = np.zeros(n_bins)
    bin_end   = np.zeros(n_bins)
    
    bin_start[0] = y_sort0[0]
    bin_end[0]   = y_sort0[r_bin_pts-1]
    for ii in range(1,n_bins):
        bin_start[ii] = bin_end[ii-1]+0.000001
        bin_end[ii]   = y_sort0[np.min([len(lO32_all)-1,(ii+1)*r_bin_pts-1])]


    #Plotting
    ctype = ['red','magenta','green','cyan','blue','black']

    if len(marker) == 0:
        marker = ['o'] * n_sample

    if fit == True:
        p0 = jiang18_coeffs

        OH_all = np.array([])
        lR23_all = np.array([])
        for nn in range(n_sample):
            OH_all   = np.append(OH_all, OH[nn])
            lR23_all = np.append(lR23_all, lR23[nn])

        opt, cov = curve_fit(O32_OH_fit, (OH_all, lO32_all), lR23_all, p0=p0)
        print(opt)

    for nn in range(n_sample):
        if len(label) != 0:
            x1 = xra[0] + 0.025*(xra[1]-xra[0])
            y1 = yra[1] - (nn*0.035 + 0.05)*(yra[1]-yra[0])
            x2 = xra[0] + 0.035*(xra[1]-xra[0])
            y2 = yra[1] - (nn*0.035 + 0.0525)*(yra[1]-yra[0])
            ax.text(x2, y2, label[nn], fontsize=8, va='center', ha='left')
            ax.plot([x1],[y1], marker=marker[nn], color='black')

        for ii in range(n_bins):
            y_ii_min = bin_start[ii] #bin_y_min + ii * dy
            y_ii_max = bin_end[ii]   #y_min + (ii+1) * dy
            idx = np.where((lO32[nn] >= y_ii_min) & (lO32[nn] <= y_ii_max))[0]

            ii_label = ''
            if nn == 0: #n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) & (lO32_all <= y_ii_max))[0]
                ii_label = r' %.2f < $\log(O_{32})$ < %.2f, N = %i' % (y_ii_min, y_ii_max,
                                                                       len(idx_all))
            if len(idx) > 0:
                ax.scatter(lR23[nn][idx], OH[nn][idx], color=ctype[ii], marker=marker[nn],
                           alpha=0.5, label=ii_label)

            if len(OH_err) != 0:
                ax.errorbar(lR23[nn][idx], OH[nn][idx], yerr=OH_err[nn][:,idx],
                            mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5,
                            fmt=None, label=None)

            if len(lR23_err) != 0:
                ax.errorbar(lR23[nn][idx], OH[nn][idx], xerr=lR23_err[nn][:,idx],
                            mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5,
                            fmt=None, label=None)

            if nn == 0: #n_sample-1:
                lO32_avg = np.average(lO32_all[idx_all])
                if fit == False:
                    opt = jiang18_coeffs

                mod_logR23 = O32_OH_fit((x_arr, lO32_avg), *opt)
                ax.annotate('%.2f' % lO32_avg, [mod_logR23[-1], x_arr[-1]],
                            color=ctype[ii], xycoords='data', ha='center',
                            va='bottom', fontsize=8)

                ax.plot(mod_logR23, x_arr, color=ctype[ii], linestyle='dashed')
        #endfor
    #endfor

    if len(xra) != 0: ax.set_xlim(xra)
    if len(yra) != 0: ax.set_ylim(yra)

    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})$')
    leg = ax.legend(loc='lower right', scatterpoints=1, fontsize=8, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.075, right=0.99, bottom=0.08, top=0.97)
    fig.savefig(out_pdf)

    # Plot differences between model and data
    if fit == False:
        out_diff_pdf = out_pdf.replace('.pdf', '.diff.pdf')
        plot_differences(lR23, lO32, OH, lO32_all, out_diff_pdf, bin_start, bin_end, n_bins=n_bins,
                         lR23_err=lR23_err, OH_err=OH_err, OH_range=yra,
                         dR23_range=dR23_range, marker=marker, label=label)
        if silent == False: log.info('### End main : '+systime())
#enddef

def get_measurements(data):

    lR23 = np.log10(data['R23'].data)
    lO32 = np.log10(data['O32'].data)
    OH   = data['OH'].data

    OH_err   = np.row_stack((data['OH_lo'].data,data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo,lR23_hi))

    return lR23, lO32, OH, OH_err, lR23_err
#enddef

def get_zcalbase_sample(prefix, dir=''):
    dir0 = '/Users/cly/data/Metallicity/Others/Te_Repository/'
    if dir == '':
        flux_file = '%s%s/%s_sample.det4363.int.fits' % (dir0, prefix, prefix)
        Te_file = '%s%s/%s_Te_table.fits' % (dir0, prefix, prefix)
    else:
        flux_file = '%s%s/%s_sample.det4363.int.fits' % (dir0, dir, prefix)
        Te_file = '%s%s/%s_Te_table.fits' % (dir0, dir, prefix)

    if 'SDSS' in prefix:
        path_SDSS = '/Users/cly/data/Metallicity/Others/SDSS/'
        flux_file = path_SDSS + 'SDSS_DR7_det4363.OIIfix.int.fits'
        Te_file   = path_SDSS + 'SDSS_DR7_det4363_Te_table.dered.fits'

    data0 = fits.getdata(flux_file)
    data1 = fits.getdata(Te_file)

    lR23 = np.log10((data0.OII_3727_FLUX + data0.OIII_5007_FLUX) / data0.H_BETA_FLUX)
    lO32 = np.log10(data0.OIII_5007_FLUX / data0.OII_3727_FLUX)

    OH = data1.OH_gas
    OH_err = np.row_stack((data1.OH_gas_lo,data1.OH_gas_hi))

    return lR23, lO32, OH
#enddef

def get_DEEP2(path0):
    infile = path0 + 'DEEP2_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err = get_measurements(data)

    return data, lR23, lO32, OH, OH_err, lR23_err
#enddef

def get_MACT(path0):
    infile = path0 + 'MACT_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    dup = ['Keck10', 'Keck17', 'Keck22', 'Keck25']
    d_idx1, d_idx2 = match_nosort_str(data['ID'].data, dup)

    data.remove_rows(d_idx1)

    lR23, lO32, OH, OH_err, lR23_err = get_measurements(data)

    return data, lR23, lO32, OH, OH_err, lR23_err
#enddef

def DEEP2_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err = get_DEEP2(path0)

    out_pdf = path0 + 'DEEP2_R23_O32_Jiang18.pdf'
    main([lR23], [lO32], [OH], out_pdf, lR23_err=[lR23_err], OH_err=[OH_err],
         xra=[0.75,1.05], yra=[7.1,8.65])

    # Got RuntimeError
    #out_pdf = path0 + 'DEEP2_R23_O32_Jiang18.fit.pdf'
    #main([lR23], [lO32], [OH], out_pdf, lR23_err=[lR23_err], OH_err=[OH_err],
    #     xra=[0.75,1.05], yra=[7.1,8.65], fit=True)

#enddef

def MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err = get_MACT(path0)

    out_pdf = path0 + 'MACT_R23_O32_Jiang18.pdf'
    main([lR23], [lO32], [OH], out_pdf, n_bins=6, lR23_err=[lR23_err],
         OH_err=[OH_err], xra=[0.60,1.15], yra=[7.10,8.7])

    out_pdf = path0 + 'MACT_R23_O32_Jiang18.fit.pdf'
    main([lR23], [lO32], [OH], out_pdf, n_bins=6, lR23_err=[lR23_err],
         OH_err=[OH_err], xra=[0.60,1.15], yra=[7.10,8.7], fit=True)

#enddef

def DEEP2_MACT_OIII4363(include_stack=False, fit=False):

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    # DEEP2
    DEEP2_data, DEEP2_lR23, DEEP2_lO32, DEEP2_OH, DEEP2_OH_err, \
        DEEP2_lR23_err = get_DEEP2(path0)

    # MACT
    MACT_data, MACT_lR23, MACT_lO32, MACT_OH, MACT_OH_err, \
        MACT_lR23_err = get_MACT(path0)

    # Combine together
    lR23 = [MACT_lR23, DEEP2_lR23]
    lO32 = [MACT_lO32, DEEP2_lO32]
    OH   = [MACT_OH,   DEEP2_OH]

    OH_err   = [MACT_OH_err, DEEP2_OH_err]
    lR23_err = [MACT_lR23_err, DEEP2_lR23_err]

    # Include RLeimbach stacked results
    if include_stack:
        stack_infile  = '/Users/cly/Google Drive/Zcalbase_gal/' + \
                        'Double_Bin_temperatures_metalicity.tbl'
        log.info('### Reading : '+stack_infile)
        stack_tbl = asc.read(stack_infile)

        stack_det = np.where((stack_tbl['S/N_4363'] >= 3.0) &
                             (stack_tbl['com_O_log'] >= 8.1))[0]
        print('stack_det : ', len(stack_det))
        R23_stack = stack_tbl['R23_Composite'].data[stack_det]
        O32_stack = stack_tbl['O32_Composite'].data[stack_det]
        OH_stack  = stack_tbl['com_O_log'].data[stack_det]

        lR23 += [R23_stack]
        lO32 += [O32_stack]
        OH   += [OH_stack]

        OH_err   += [np.array([[0] * len(stack_det), [0] * len(stack_det)])]
        lR23_err += [np.array([[0] * len(stack_det), [0] * len(stack_det)])]
        print(OH_err)
    #endif

    out_pdf = path0 + 'MACT_DEEP2_R23_O32_Jiang18.pdf'
    label = [r'$\mathcal{MACT}$  (Ly+2016)',
             r'DEEP2 [OIII]$\lambda$4363-detected']
    marker = ['o', '*']

    if include_stack:
        out_pdf = out_pdf.replace('.pdf', '.stack.pdf')
        label += ['DEEP2 stacked detections']
        marker += ['s']
        yra = [7.10,9.0]
    else:
        yra = [7.10,8.8]

    main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.45,1.15], yra=yra, dR23_range=[-0.15,0.3], marker=marker, label=label)

    if fit:
        out_pdf = path0 + 'MACT_DEEP2_R23_O32_Jiang18.fit.pdf'
        main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
             xra=[0.45,1.15], yra=yra, marker=['*','o'], label=label, fit=True)

def Zcalbase():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    ref_name0 = ['Berg2012','Kennicutt2003','Izotov1994','Thuan1995',
                 'Izotov1997','Guseva2009', 'Izotov2012', 'SDSS']
    dir0      = ['','','BCGs','BCGs','Pilyugin2012/Izotov1997',
                 'Pilyugin2012/Guseva2009', 'Pilyugin2012/Izotov2012','']

    lR23_all = []
    lO32_all = []
    OH_all   = []

    for name,dir in zip(ref_name0,dir0):
        lR23, lO32, OH = get_zcalbase_sample(name, dir=dir)

        lR23_all.append(lR23)
        lO32_all.append(lO32)
        OH_all.append(OH)

    out_pdf = path0 + 'Zcalbase_Jiang18.pdf'
    label = ref_name0 #['Kennicutt+2003']
    main(lR23_all, lO32_all, OH_all, out_pdf, n_bins=6, xra=[0.1,1.05],
         yra=[7.0,9.0], marker=['s','x','s','s','D','D','s','*'], label=label)
