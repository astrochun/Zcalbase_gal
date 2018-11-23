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

def main(lR23, lO32, OH, out_pdf, silent=False, verbose=True):

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

    y_min = np.min(lO32)
    y_max = np.max(lO32)

    ax.scatter(lR23, OH)

    if silent == False: log.info('### End main : '+systime())
#enddef

