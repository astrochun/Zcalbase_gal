"""
check_astrometry
================

Compares the WCS for two different catalogs, such as a target sample and a
reference survey (e.g., 2MASS, SDSS) to check relative astrometry
"""

import sys, os

from chun_codes import systime
from chun_codes import match_nosort
from chun_codes import chun_crossmatch

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
    
def main(infile1, infile2=None, cat2_survey='SDSS', out_pdf=None, silent=False,
         verbose=True):

    '''
    Main function to cross-match two catalogs and plot relative astrometric
    accuracy. The second input catalog is consider the reference catalog
    (i.e., the one with accurate astrometry)

    Parameters
    ----------
    infile1 : string
      Filename for input catalog file

    infile2 : string
      Filename for input reference catalog file. Default: None
      Either provide this or cat2_survey to query appropriate servers

    cat2_survey : string
      Name of survey to query. Default: SDSS
      Either provide this keyword or infile2 but not both

    out_pdf : string
      Output PDF file name. Default: None -> based on infile1 name
      
    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 6 January 2017
    '''
    
    if silent == False: print '### Begin check_astrometry.main | '+systime()

    if out_pdf == None:
        out_pdf = infile1.replace('.fits.gz','.pdf').replace('.dat','.pdf')
    if silent == False: print '### End check_astrometry.main | '+systime()
#enddef

