"""
locate_em_lines
===============

Generate plots illustrating location of nebular emission lines given redshift,
the OH night skyline spectrum, and the atmospheric transmission
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

def main(in_cat, out_pdf, lambda0_min, lambda0_max, format='commented_header',
         silent=False, verbose=True):

    '''
    Main function to read in catalog of sources and output PDF file
    illustrating emission line location, OH night skyline spectrum,
    and atmospheric transmission

    Parameters
    ----------
    in_cat : string
      Filename for ASCII catalog containing source name and redshift

    out_pdf : string
      Filename to output PDF. Full path should be provided

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 17 December 2016
    '''

    if silent == False: print '### Begin locate_em_lines.main() | '+systime()

    if silent == False: print '### Reading : ', in_cat
    cat_data = asc.read(in_cat)

    ID    = cat_data['ID']
    zspec = cat_data['redshift']

    # Values are in Angstrom
    OII_loc   = (1.0+zspec) * 3727.00
    Hb_loc    = (1.0+zspec) * 4861.32
    OIIIa_loc = (1.0+zspec) * 4958.91
    OIIIb_loc = (1.0+zspec) * 5006.84
    Ha_loc    = (1.0+zspec) * 6562.80
    NIIa_loc  = (1.0+zspec) * 6548.10
    NIIb_loc  = (1.0+zspec) * 6583.60
    SIIa_loc  = (1.0+zspec) * 6717.42
    SIIb_loc  = (1.0+zspec) * 6730.78

    if silent == False: print '### End locate_em_lines.main() | '+systime()
#enddef

