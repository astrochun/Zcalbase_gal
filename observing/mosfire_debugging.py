"""
mosfire_debugging
=================

Some troubleshooting for when MOSFIRE DRP breaks
"""

import sys, os

from chun_codes import systime

from os.path import exists

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength, Extract

flatops = Options.flat
waveops = Options.wavelength

def main(path, maskname, band, silent=False, verbose=True):

    '''
    Main function for testing

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
    Created by Chun Ly, 18 June 2018
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    
    if silent == False: log.info('### End main : '+systime())
#enddef

