"""
sdss_2mass_proper_motion
========================

Function to determine proper motion based on 2MASS and SDSS coordinates

Requirements:
 astroquery
 matplotlib
"""

import sys, os

from chun_codes import systime

from os.path import exists

import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u

from astropy.table import Table, Column

from astropy import coordinates as coords
from astroquery.sdss import SDSS
from astroquery.irsa import Irsa as IRSA

from astropy.time import Time

from scipy.stats import linregress

SDSS_fld      = ['ra','dec','raErr','decErr','objid','run','rerun','camcol',
                 'field','obj', 'type','mode','mjd']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def main(tab0, silent=True, verbose=False):
    '''
    Main function to determine proper motion based on 2MASS and SDSS coordinates

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Table containing alignment stars SDSS and 2MASS information

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------
    pra0, pdec0 : float or array like
      Proper motion in mas/yr for RA and Dec

    Notes
    -----
    Created by Chun Ly, 23 January 2017
     - Started as a copy of observing.find_nearby_bright_star.sdss_2mass_proper_motion
     - Decided to make it a standalone and to run on full sample
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.main() | '+systime()

    len0 = len(tab0)

    with_2mass = np.array([idx for idx in range(len(tab0)) if
                           'XXXX' not in tab0['date_2mass'][idx]])

    # Initialize | + on 22/01/2017
    pra0  = np.repeat(-999.999, len0)
    pdec0 = np.repeat(-999.999, len0)

    # Deal with those without 2MASS data | + on 22/01/2017
    if len(with_2mass) > 0:
        tab2 = tab0.copy()
        tab2 = tab2[with_2mass]
        c_sdss  = coords.SkyCoord(ra=tab2['ra'], dec=tab2['dec'], unit=(u.deg))
        c_2mass = coords.SkyCoord(ra=tab2['ra_2mass'], dec=tab2['dec_2mass'],
                                  unit=(u.deg))

        ra_2mass  = c_2mass.ra.value
        dec_2mass = c_2mass.dec.value

        for ii in range(len(tab2)):
            tab_SDSS = SDSS.query_region(c_sdss[ii], radius=2*u.arcsec,
                                         data_release=12,
                                         photoobj_fields=SDSS_phot_fld)
            SDSS_date = Time(tab_SDSS['mjd'].data, format='mjd')
            SDSS_date.format='decimalyear'

            SDSS_ra     = tab_SDSS['ra'].data
            SDSS_raErr  = tab_SDSS['raErr'].data
            SDSS_dec    = tab_SDSS['dec'].data
            SDSS_decErr = tab_SDSS['decErr'].data

            t_2mass = Time(tab2['date_2mass'][ii], format='iso')
            t_2mass.format='decimalyear'

            time0 = np.array([t_2mass.value] + list(SDSS_date.value))
            ra0   = [ra_2mass[ii]]  + SDSS_ra.tolist()
            dec0  = [dec_2mass[ii]] + SDSS_dec.tolist()

            x0 = time0 - 2000.0
            ra_regress  = linregress(x0, ra0)
            dec_regress = linregress(x0, dec0)
            pra0[with_2mass[ii]]  = ra_regress.slope * 3600.0 * 1E3
            pdec0[with_2mass[ii]] = dec_regress.slope * 3600.0 * 1E3
        #endfor
    #endif

    if silent == False:
        print '### End sdss_2mass_proper_motion.main() | '+systime()

    return pra0, pdec0
#enddef

def old(tab0, silent=True, verbose=False):
    '''
    Function to determine proper motion based on 2MASS and SDSS coordinates

    Parameters
    ----------
    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------
    pra0, pdec0 : float or array like
      Proper motion in mas/yr for RA and Dec

    Notes
    -----
    Created by Chun Ly, 21 January 2017
    Modified by Chun Ly, 22 January 2017
     - Fix small bug
     - Handle cases without 2MASS data
     - Small bug found. Require [with_2mass] as non empty
    This code came from Zcalbase_gal.observing.find_nearby_bright_star.
    It is now obsolete (replaced with main() function above)
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.old | '+systime()

    with_2mass = np.where('XXXX' not in tab0['date_2mass'])[0] # + on 22/01/2017

    # Initialize | + on 22/01/2017
    pra0  = np.repeat(0.000, len(tab0))
    pdec0 = np.repeat(0.000, len(tab0))

    # Deal with those without 2MASS data | + on 22/01/2017
    if len(with_2mass) > 0:
        tab2 = tab0.copy()
        tab2 = tab2[with_2mass]
        c_sdss  = coords.SkyCoord(ra=tab2['ra'], dec=tab2['dec'], unit=(u.deg))
        c_2mass = coords.SkyCoord(ra=tab2['ra_2mass'], dec=tab2['dec_2mass'],
                                  unit=(u.deg))

        dra0, ddec0 = c_2mass.spherical_offsets_to(c_sdss)

        t_sdss  = Time(tab2['mjd'], format='mjd')
        t_2mass = Time(tab2['date_2mass'], format='iso') # Mistake fix on 22/01/2017
        t_diff  = (t_sdss - t_2mass).to(u.yr).value # in years

        # Mod on 22/01/2017
        pra0[with_2mass]  = dra0.to(u.mas).value/t_diff
        pdec0[with_2mass] = ddec0.to(u.mas).value/t_diff
    #endif

    if silent == False:
        print '### End sdss_2mass_proper_motion.old | '+systime()

    return pra0, pdec0
#enddef
