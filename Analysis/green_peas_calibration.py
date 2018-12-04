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

from astropy.table import Table
from astropy import log

jiang18_coeffs = [-24.135, 6.1523, -0.37866, -0.147, -7.071]

def O32_OH_fit(x, y, a, b, c, d, e):
    '''
    Parameters
    ----------
     x : 12+log(O/H)
     y : log([OIII]/[OII])
    '''
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

    logR23 = O32_OH_fit(x, y, *jiang18_coeffs)

    return logR23
#enddef

def main(lR23, lO32, OH, out_pdf, n_bins=4, lR23_err=[], OH_err=[], xra=[],
         yra=[], marker=[], label=[], fit=False, silent=False, verbose=True):

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

        OH_min1[nn] = np.min(OH[nn])
        OH_max1[nn] = np.max(OH[nn])

    # Grid of 12+log(O/H)
    x_arr = np.arange(min(OH_min1),max(OH_max1),0.05)

    y_min = np.min(min1)
    y_max = np.max(max1)

    lO32_all = np.array([])
    for nn in range(n_sample):
        lO32_all = np.append(lO32_all, lO32[nn])

    sort0   = np.argsort(lO32_all)
    y_sort0 = lO32_all[sort0]

    bin_pts = np.int(len(lO32_all)/n_bins)

    bin_start = np.zeros(n_bins)
    bin_end   = np.zeros(n_bins)
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*bin_pts]
        bin_end[ii]   = y_sort0[(ii+1)*bin_pts-1]

    ctype = ['red','magenta','green','cyan','blue','black']

    if len(marker) == 0:
        marker = ['o'] * n_sample

    if fit == True:
        p0 = jiang18_coeffs
        opt, cov = curve_fit(O32_OH_fit, OH, lO32, lR23, p0=p0)

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
            if nn == n_sample-1:
                idx_all = np.where((lO32_all >= y_ii_min) & (lO32_all <= y_ii_max))[0]
                ii_label = r' %.2f < $\log(O_{32})$ < %.2f, N = %i' % (y_ii_min, y_ii_max,
                                                                       len(idx_all))
            if len(idx) > 0:
                ax.scatter(lR23[nn][idx], OH[nn][idx], color=ctype[ii], marker=marker[nn],
                           alpha=0.5, label=ii_label)

            if len(OH_err[nn]) != 0:
                ax.errorbar(lR23[nn][idx], OH[nn][idx], yerr=OH_err[nn][:,idx],
                            mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5,
                            fmt=None, label=None)

            if len(lR23_err[nn]) != 0:
                ax.errorbar(lR23[nn][idx], OH[nn][idx], xerr=lR23_err[nn][:,idx],
                            mec=ctype[ii], ecolor=ctype[ii], capsize=0, alpha=0.5,
                            fmt=None, label=None)

            if nn == n_sample-1:
                lO32_avg = np.average(lO32_all[idx_all])
                j18_logR23 = jiang18(x_arr, lO32_avg)
                ax.annotate('%.2f' % lO32_avg, [j18_logR23[-1], x_arr[-1]],
                            color=ctype[ii], xycoords='data', ha='center',
                            va='bottom', fontsize=8)

                ax.plot(j18_logR23, x_arr, color=ctype[ii], linestyle='dashed')
        #endfor
    #endfor

    if len(xra) != 0: ax.set_xlim(xra)
    if len(yra) != 0: ax.set_ylim(yra)

    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})$')
    leg = ax.legend(loc='lower right', scatterpoints=1, fontsize=8, framealpha=0.5)
    for lh in leg.legendHandles:
        lh.set_alpha(0.5)

    plt.subplots_adjust(left=0.075, right=0.99, bottom=0.08, top=0.99)
    fig.savefig(out_pdf)

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

#enddef

def MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err = get_MACT(path0)

    out_pdf = path0 + 'MACT_R23_O32_Jiang18.pdf'
    main([lR23], [lO32], [OH], out_pdf, n_bins=6, lR23_err=[lR23_err],
         OH_err=[OH_err], xra=[0.60,1.15], yra=[7.10,8.7])

#enddef

def DEEP2_MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    # DEEP2
    DEEP2_data, DEEP2_lR23, DEEP2_lO32, DEEP2_OH, DEEP2_OH_err, \
        DEEP2_lR23_err = get_DEEP2(path0)

    # MACT
    MACT_data, MACT_lR23, MACT_lO32, MACT_OH, MACT_OH_err, \
        MACT_lR23_err = get_MACT(path0)

    # Combine together
    lR23 = [DEEP2_lR23, MACT_lR23]
    lO32 = [DEEP2_lO32, MACT_lO32]
    OH   = [DEEP2_OH,   MACT_OH]

    OH_err = [DEEP2_OH_err, MACT_OH_err]

    lR23_err = [DEEP2_lR23_err, MACT_lR23_err]

    out_pdf = path0 + 'MACT_DEEP2_R23_O32_Jiang18.pdf'
    label = [r'DEEP2 [OIII]$\lambda$4363-detected',
             r'$\mathcal{MACT}$  (Ly+2016)']
    main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.6,1.15], yra=[7.10,8.7], marker=['*','o'], label=label)

