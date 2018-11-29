"""
local_analog_calibration
====

Plots R23 and O32 line ratios and compare it to Bian et al. (2018)
calibration that is based on local analogs to z~2 galaxies
"""

import sys, os

from chun_codes import systime, match_nosort_str

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

def bian18_R23(OH):
    '''
    Function to return log(R23) given metallicity

    Parameters
    ----------
     OH : 12+log(O/H)
    '''

    R23_coeff = [138.0430, -54.8284, 7.2954, -0.32293]
    R23_p = np.poly1d(R23_coeff)

    return R23_p(OH)
#enddef

def bian18_O32(O32):
    '''
    Function to return metallicity given log(O32)

    Parameters
    ----------
     O32 : log([OIII]/[OII])
    '''

    OH = 8.54 - 0.59 * O32

    return OH
#enddef

def main(lR23, lO32, OH, out_pdf, n_bins=4, lR23_err=[], OH_err=[],
         R23_xra=[], yra=[], silent=False, verbose=True):

    '''
    Main function to plot dataset against Bian+ (2018) calibration

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
    Created by Chun Ly, 29 November 2018
    '''

    if silent == False: log.info('### Begin main : '+systime())

    fig, ax = plt.subplots(n_cols=2)

    O32_min = np.min(lO32)
    O32_max = np.max(lO32)
    O32_arr = np.arange(O32_min,O32_max, 0.025)

    bian_OH = bian18_O32(O32_arr)

    # Grid of 12+log(O/H)
    OH_arr = np.arange(min(OH),max(OH),0.05)

    bian_R23 = bian18_R23(OH_arr)
    
    ax[0].scatter(lR23, OH, color='blue', marker='o', alpha=0.5)
    ax[0].plot(bian_R23, OH_arr)

    if len(OH_err) != 0:
        ax[0].errorbar(lR23, OH, yerr=OH_err, mec='blue', ecolor='blue',
                       capsize=0, alpha=0.5, fmt=None, label=None)

    if len(lR23_err) != 0:
        ax[0].errorbar(lR23, OH, xerr=lR23_err, mec='blue', ecolor='blue',
                       capsize=0, alpha=0.5, fmt=None, label=None)

    if len(R23_xra) != 0: ax[0].set_xlim(R23_xra)
    if len(yra) != 0: ax[0].set_ylim(yra)

    ax[0].set_xlabel(r'$\log(R_{23})$')
    ax[0].set_ylabel(r'$12+\log({\rm O/H})$')

    ax[1].scatter(lO32, OH, color='blue', marker='o', alpha=0.5)
    if len(OH_err) != 0:
        ax[0].errorbar(lO32, OH, yerr=OH_err, mec='blue', ecolor='blue',
                       capsize=0, alpha=0.5, fmt=None, label=None)

    if len(lO32_err) != 0:
        ax[0].errorbar(lO32, OH, xerr=lR23_err, mec='blue', ecolor='blue',
                       capsize=0, alpha=0.5, fmt=None, label=None)


    plt.subplots_adjust(left=0.075, right=0.99, bottom=0.08, top=0.99)
    fig.savefig(out_pdf)

    if silent == False: log.info('### End main : '+systime())
#enddef
