"""
green_peas_calibration
====

Plots R23 and O32 line ratios and compare it to Jiang et al. (2018)
calibration that is based on green pea galaxies
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

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

def main(silent=False, verbose=True):

    '''
    Main function for green_peas_calibration

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, XX
    '''

    
    if silent == False: log.info('### Begin main : '+systime())
    
    if silent == False: log.info('### End main : '+systime())
#enddef

