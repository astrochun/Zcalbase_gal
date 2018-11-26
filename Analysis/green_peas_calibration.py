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

from astropy.table import Table
from astropy import log

def jiang18(x, y):
    '''
    Function to return log(R23) based on metallicity and [OIII]/[OII] flux ratio

    Parameters
    ----------
     x : 12+log(O/H)
     y : log([OIII]/[OII])
    '''

    a = -24.135
    b =   6.1523
    c =  -0.37866
    d =  -0.147
    e =  -7.071

    logR23 = a + b*x + c * x**2 - d * (e + x) * y

    return logR23
#enddef

def main(lR23, lO32, OH, out_pdf, n_bins=4, lR23_err=[], OH_err=[],
         xra=[], yra=[], silent=False, verbose=True):

    '''
    Main function to plot dataset against Jiang+ (2018) calibration

    Parameters
    ----------
    lR23 : log(R_23)
    lO32 : log([OIII]/[OII])
    OH   : 12+log(O/H)

    out_pdf : full path for output PDF

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

    # Grid of 12+log(O/H)
    x_arr = np.arange(min(OH),max(OH),0.05)

    y_min = np.min(lO32)
    y_max = np.max(lO32)

    sort0   = np.argsort(lO32)
    y_sort0 = lO32[sort0]

    bin_pts = np.int(len(lO32)/n_bins)

    bin_start = np.zeros(n_bins)
    bin_end   = np.zeros(n_bins)
    for ii in range(n_bins):
        bin_start[ii] = y_sort0[ii*bin_pts]
        bin_end[ii]   = y_sort0[(ii+1)*bin_pts-1]

    ctype = ['red','magenta','green','cyan','blue','black']

    for ii in range(n_bins):
        y_ii_min = bin_start[ii] #bin_y_min + ii * dy
        y_ii_max = bin_end[ii]   #y_min + (ii+1) * dy
        idx = np.where((lO32 >= y_ii_min) & (lO32 <= y_ii_max))[0]
        ii_label = r' %.2f < $\log(O_{32})$ < %.2f, N = %i' % (y_ii_min, y_ii_max,
                                                               len(idx))
        if len(idx) > 0:
            ax.scatter(lR23[idx], OH[idx], color=ctype[ii], marker='o', alpha=0.5,
                       label=ii_label)

            if len(OH_err) != 0:
                ax.errorbar(lR23[idx], OH[idx], yerr=OH_err[:,idx], mec=ctype[ii],
                            ecolor=ctype[ii], capsize=0, alpha=0.5, fmt=None, label=None)

            if len(lR23_err) != 0:
                ax.errorbar(lR23[idx], OH[idx], xerr=lR23_err[:,idx], mec=ctype[ii],
                            ecolor=ctype[ii], capsize=0, alpha=0.5, fmt=None, label=None)

            lO32_avg = np.average(lO32[idx])
            j18_logR23 = jiang18(x_arr, lO32_avg)
            ax.annotate('%.2f' % lO32_avg, [j18_logR23[-1], x_arr[-1]],
                        color=ctype[ii], xycoords='data', ha='center',
                        va='bottom', fontsize=8)

            ax.plot(j18_logR23, x_arr, color=ctype[ii], linestyle='dashed')

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

def DEEP2_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    infile = path0 + 'DEEP2_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    lR23 = np.log10(data['R23'])
    lO32 = np.log10(data['O32'])
    OH   = data['OH']

    OH_err   = np.row_stack((data['OH_lo'].data,data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo,lR23_hi))

    #print np.min(lR23), np.max(lR23)
    #print np.min(OH), np.max(OH)
    out_pdf = path0 + 'DEEP2_R23_O32_Jiang18.pdf'
    main(lR23, lO32, OH, out_pdf, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.75,1.05], yra=[7.1,8.65])

#enddef

def MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    infile = path0 + 'MACT_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    dup = ['Keck10', 'Keck17', 'Keck22', 'Keck25']
    d_idx1, d_idx2 = match_nosort_str(data['ID'].data, dup)

    data.remove_rows(d_idx1)

    lR23 = np.log10(data['R23'])
    lO32 = np.log10(data['O32'])
    OH   = data['OH']

    OH_err = np.row_stack((data['OH_lo'].data,data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo,lR23_hi))

    out_pdf = path0 + 'MACT_R23_O32_Jiang18.pdf'
    main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.60,1.15], yra=[7.10,8.7])

#enddef

def DEEP2_MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    # DEEP2

    infile = path0 + 'DEEP2_R23_O32_derived.tbl'
    log.info('### Reading : '+infile)
    deep2_data = asc.read(infile)

    deep2_lR23 = np.log10(deep2_data['R23'])
    deep2_lO32 = np.log10(deep2_data['O32'])
    deep2_OH   = deep2_data['OH']

    deep2_lR23_lo = deep2_lR23 - np.log10(deep2_data['R23'] - deep2_data['R23_lo'])
    deep2_lR23_hi = np.log10(deep2_data['R23'] + deep2_data['R23_hi']) - deep2_lR23

    # MACT
    
    infile = path0 + 'MACT_R23_O32_derived.tbl'
    log.info('### Reading : '+infile)
    mact_data = asc.read(infile)

    dup = ['Keck10', 'Keck17', 'Keck22', 'Keck25']
    d_idx1, d_idx2 = match_nosort_str(mact_data['ID'].data, dup)
    mact_data.remove_rows(d_idx1)

    mact_lR23 = np.log10(mact_data['R23'])
    mact_lO32 = np.log10(mact_data['O32'])
    mact_OH   = mact_data['OH']

    mact_lR23_lo = mact_lR23 - np.log10(mact_data['R23'] - mact_data['R23_lo'])
    mact_lR23_hi = np.log10(mact_data['R23'] + mact_data['R23_hi']) - mact_lR23

    # Combine together
    lR23 = np.concatenate((deep2_lR23, mact_lR23))
    lO32 = np.concatenate((deep2_lO32, mact_lO32))
    OH   = np.concatenate((deep2_OH,   mact_OH))

    OH_lo = np.concatenate((deep2_data['OH_lo'].data, mact_data['OH_lo'].data))
    OH_hi = np.concatenate((deep2_data['OH_hi'].data, mact_data['OH_hi'].data))
    OH_err = np.row_stack((OH_lo,OH_hi))

    lR23_lo = np.concatenate((deep2_lR23_lo, mact_lR23_lo))
    lR23_hi = np.concatenate((deep2_lR23_hi, mact_lR23_hi))
    lR23_err = np.row_stack((lR23_lo, lR23_hi))

    out_pdf = path0 + 'MACT_DEEP2_R23_O32_Jiang18.pdf'
    main(lR23, lO32, OH, out_pdf, n_bins=6, lR23_err=lR23_err, OH_err=OH_err,
         xra=[0.6,1.15], yra=[7.10,8.7])

