"""
check_astrometry
================

Compares the WCS for two different catalogs, such as a target sample and a
reference survey (e.g., 2MASS, SDSS) to check relative astrometry
"""

import sys, os

from chun_codes import systime
from chun_codes import chun_crossmatch

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table

import astropy.coordinates as coords
import astropy.units as u

def main(infile1, c0, coords=['ALPHA_J2000', 'DELTA_J2000'], infile2=None,
         cat2_survey='SDSS', out_pdf=None, silent=False, verbose=True):

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
    Modified by Chun Ly, 12 January 2017
     - Read in infile1 based on format type (ASCII or FITS)
    '''
    
    if silent == False: print '### Begin check_astrometry.main | '+systime()

    # + on 12/01/2017
    if '.fits' in infile1:
        hdu1  = fits.open(infile1)
        data1 = hdu1[1].data
    else:
        data1 = asc.read(infile1)

    if out_pdf == None:
        out_pdf = infile1.replace('.fits.gz','.pdf').replace('.dat','.pdf')
    if silent == False: print '### End check_astrometry.main | '+systime()
#enddef

def SDF(silent=False, verbose=True):
    '''
    Cross-match SDF NB catalogs against SDSS to get astrometry around each
    [OIII]4363 targets from Ly et al. (2016), ApJS, 226, 5

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
    Created by Chun Ly, 12 January 2017
    '''
    
    if silent == False: print '### Begin check_astrometry.SDF | '+systime()

    path0     = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    cat_path0 = '/Users/cly/data/SDF/NBcat/SEx/'

    infile0 = path0 + 'targets_sdf_astro.2017a.txt'
    if silent == False: print '### Reading : ', infile0
    data0   = asc.read(infile0, format='commented_header')
    print data0
    
    n_sources = len(data0)

    ID  = data0['ID']
    RA  = data0['RA']
    DEC = data0['DEC']
    c0  = coords.SkyCoord(ra=RA, dec=DEC, unit=(u.hour, u.deg))
    
    cat_files = [cat_path0+f0+'emitters_allphot.fits' for f0 in data0['filt']]
    print cat_files

    box_size = (5 * u.arcmin).to(u.arcsec).value
    for ii in range(n_sources):
        if silent == False: print '### Reading : ', cat_files[ii]
        data1 = fits.getdata(cat_files[ii])
        RA1   = data1.ALPHA_J2000
        DEC1  = data1.DELTA_J2000

        t_c = c0[ii]
        RA_off  = (RA1 - t_c.ra.deg) * 3600.0 * np.cos(np.radians(t_c.dec.deg))
        DEC_off = (DEC1 - t_c.dec.deg) * 3600.0
        in_box = np.where((np.absolute(RA_off) <= box_size) &
                          (np.absolute(DEC_off) <= box_size))[0]
        box_data1 = data1[in_box]
        
        if silent == False: print '## in_box : ', len(in_box)
    #endfor
    if silent == False: print '### End check_astrometry.SDF | '+systime()
#enddef
