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
from matplotlib.backends.backend_pdf import PdfPages

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

    # Values are in Angstroms
    co_filename = __file__
    lambda0_file = os.path.dirname(co_filename) + '/' + 'lambda0.txt'
    lambda0_data = asc.read(lambda0_file, format='commented_header')
    
    lambda0 = lambda0_data['lambda0']
    #[3727, 4861.32, 4958.91, 5006.84, 6562.80, 6548.10, 6583.60, 6717.42, 6730.78]
    
    #OII_loc   = (1.0+zspec) * 3727.00
    #Hb_loc    = (1.0+zspec) * 4861.32
    #OIIIa_loc = (1.0+zspec) * 4958.91
    #OIIIb_loc = (1.0+zspec) * 5006.84
    #Ha_loc    = (1.0+zspec) * 6562.80
    #NIIa_loc  = (1.0+zspec) * 6548.10
    #NIIb_loc  = (1.0+zspec) * 6583.60
    #SIIa_loc  = (1.0+zspec) * 6717.42
    #SIIb_loc  = (1.0+zspec) * 6730.78

    n_sources = len(ID)
    
    n_panels = 4 # Number of figure panels

    pp = PdfPages(out_pdf)

    for ii in range(n_sources):
        c_ii = ii % n_panels
        if c_ii == 1: fig, (ax1, ax2, ax3, ax4) = plt.subplots(4) #, sharex=True)

        cmd1 = 'ax0 = ax'+str(c_ii)
        execute(cmd1)

        ax0.plot()
        ax0.plot(Ha_loc[ii]
        
    if silent == False: print '### End locate_em_lines.main() | '+systime()
#enddef

