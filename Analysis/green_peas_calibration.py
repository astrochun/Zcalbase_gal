"""
green_peas_calibration
====

Plots R23 and O32 line ratios and compare it to Jiang et al. (2018)
calibration that is based on green pea galaxies
"""

import sys, os

from chun_codes import systime

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

def main(lR23, lO32, OH, out_pdf, n_bins=5, silent=False, verbose=True):

    '''
    Main function to plot dataset against Jiang+ (2018) calibration

    Parameters
    ----------
    lR23 : log(R_23)
    lO32 : log([OIII]/[OII])
    OH   : 12+log(O/H)

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

    dy = (y_max-y_min)/n_bins

    ctype = ['red','magenta','green','cyan','blue']

    for ii in range(n_bins):
        y_ii_min = y_min + ii * dy
        y_ii_max = y_min + (ii+1) * dy
        idx = np.where((lO32 >= y_ii_min) & (lO32 <= y_ii_max))[0]
        ii_label = r' %.2f < $\log(O_{32})$ < %.2f, N = %i' % (y_ii_min, y_ii_max,
                                                               len(idx))
        if len(idx) > 0:
            ax.scatter(lR23[idx], OH[idx], color=ctype[ii], marker='o', alpha=0.5,
                       label=ii_label)

            j18_logR23 = jiang18(x_arr, np.average(lO32[idx]))

            ax.plot(j18_logR23, x_arr, color=ctype[ii], linestyle='dashed')
    ax.set_xlabel(r'$\log(R_{23})$')
    ax.set_ylabel(r'$12+\log({\rm O/H})$')
    ax.legend(loc='upper left', fontsize=10)

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

    out_pdf = path0 + 'DEEP2_R23_O32_Jiang18.pdf'
    main(lR23, lO32, OH, out_pdf)

#enddef
