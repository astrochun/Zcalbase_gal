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

from MOSFIRE import Options, Wavelength, IO

#from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, Rectify, Wavelength, Extract

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

    flat_file = path+'combflat_2d_'+band+'.fits'

    flatlist = path+'Flat.txt'
    flatlist = IO.list_file_to_strings(flatlist)

    bpos = np.ones(92) * -1

    for fname in flatlist:
        hdr, dat, bs = IO.readmosfits(fname, flatops)
        try: bs0
        except: bs0 = bs

        if np.any(bs0.pos != bs.pos):
            print "bs0: "+str(bs0.pos)+" bs: "+str(bs.pos)
            log.warn("Barset do not seem to match")

        if hdr["filter"] != band:
            log.warn("Filter name %s does not match header filter name "
                     "%s in file %s" % (band, hdr["filter"], fname))

        for i in xrange(len(bpos)):
            b = hdr["B{0:02d}POS".format(i+1)]
            if bpos[i] == -1:
                bpos[i] = b
            else:
                if bpos[i] != b:
                    log.warn("Bar positions are not all the same in "
                             "this set of flat files")

    bs = bs0

    print bs

    if silent == False: log.info('### End main : '+systime())
#enddef

