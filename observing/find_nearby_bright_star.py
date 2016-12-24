"""
find_nearby_bright_star
=======================

Search SDSS or 2MASS catalogs for nearest bright star as offset star
or to determine PA for longslit spectroscopy

Requirements:
 astroquery (https://astroquery.readthedocs.io/en/latest/)

"""

import sys, os

from chun_codes import systime
from chun_codes import chun_crossmatch

from os.path import exists

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

import astropy.units as u
from astropy.table import Table, Column

from astropy import coordinates as coords
from astroquery.sdss import SDSS
from astroquery.irsa import Irsa as IRSA

# For SDSS only
SDSS_fld      = ['ra', 'dec', 'objid', 'run', 'rerun', 'camcol', 'field']
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMag_g', 'modelMag_r',
                            'modelMag_i', 'modelMag_z']

def main(infile, max_radius=60*u.arcsec, mag_limit=20.0, mag_filt='modelMag_i',
         catalog='SDSS', format='commented_header', silent=False, verbose=True):

    '''
    Main function to find nearby star

    Parameters
    ----------
    infile : string
      Filename for ASCII file that contains ID, RA, and DEC
      RA and DEC must be provided in sexigesimal format.
      The columns should be referred to as 'ID', 'RA', and 'DEC'

    max_radius : float
      Maximum radius. Provide with astropy.units for arcsec,
      arcmin, or degrees. Default: 60 * u.arcsec

    mag_limit : float
      Faintest source to consider in AB mag. Default: 20.0

    SDSS : boolean
      Set to True if using SDSS catalogs. Default: True
    
    TWOMASS : boolean
      Set to True if using 2MASS catalogs. Default: False

    format : string
      Format of infile ASCII file. Default: "commented_header"

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 December 2016
    '''

    if silent == False:
        print '### Begin find_nearby_bright_star.main | '+systime()

    if silent == False: print '### Reading : ', infile
    data0 = asc.read(infile)

    ID  = data0['ID'].data
    RA  = data0['RA'].data
    DEC = data0['DEC'].data
    n_sources = len(ID)

    if silent == False:
        print '## Search criteria : '
        print '## max_radius(arcsec) : ', max_radius.to(u.arcsec).value
        print '## mag_limit : ', mag_limit
        print '## filter selection : ', mag_filt
    
    for ii in range(1): #range(n_sources):
        c = coords.SkyCoord(ra=RA[ii], dec=DEC[ii], unit=(u.hour, u.degree))
        if catalog == 'SDSS':
            xid = SDSS.query_region(c, max_radius, data_release=12,
                                    photoobj_fields=SDSS_phot_fld)

            
        good = np.where(xid[mag_filt] <= mag_limit)[0]
        xid = xid[good]
        if silent == False:
            print '## Finding nearby stars for '+ID[ii]+'. '+\
                str(len(good))+' found.'
            
        c1 = coords.SkyCoord(xid['ra'], xid['dec'],
                             unit=(u.degree, u.degree))
        sep = c.separation(c1).to(u.arcsec).value
        col1 = Column(sep, name='Dist(arcsec)')
        xid.add_column(col1)
        print xid
            
    if silent == False:
        print '### End find_nearby_bright_star.main | '+systime()
#enddef

