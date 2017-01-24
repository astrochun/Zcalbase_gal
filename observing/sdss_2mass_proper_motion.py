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

from matplotlib.backends.backend_pdf import PdfPages

SDSS_fld      = ['ra','dec','raErr','decErr','objid','run','rerun','camcol',
                 'field','obj', 'type','mode','mjd']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def main(tab0, out_pdf=None, silent=False, verbose=True):
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
     - Plotting of proper motion info added in later
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
        len2 = len(tab2) # later + on 23/01/2017
        c_sdss  = coords.SkyCoord(ra=tab2['ra'], dec=tab2['dec'], unit=(u.deg))
        c_2mass = coords.SkyCoord(ra=tab2['ra_2mass'], dec=tab2['dec_2mass'],
                                  unit=(u.deg))

        ra_2mass  = c_2mass.ra.value
        dec_2mass = c_2mass.dec.value

        # later + on 23/01/2017
        if out_pdf == None:
            out_pdf = 'sdss_2mass_proper_motion.pdf'
        pp = PdfPages(out_pdf)
        n_panels = 5

        for ii in range(len2):
            row = ii % n_panels
            if row == 0: fig, ax0 = plt.subplots(5, 2)

            tab_SDSS = SDSS.query_region(c_sdss[ii], radius=2*u.arcsec,
                                         data_release=12,
                                         photoobj_fields=SDSS_phot_fld)

            # Restrict to 1"
            c_tab = coords.SkyCoord(tab_SDSS['ra'], tab_SDSS['dec'],
                                    unit=(u.deg))

            dist0 = c_sdss[ii].separation(c_tab).to(u.arcsec).value
            in_rad = np.where(dist0 <= 1.0)[0]
            tab_SDSS = tab_SDSS[in_rad]

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

            # Mod on 23/01/2017
            ra_fit  = linregress(x0, ra0)
            dec_fit = linregress(x0, dec0)
            pra0[with_2mass[ii]]  = ra_fit.slope * 3600.0 * 1E3
            pdec0[with_2mass[ii]] = dec_fit.slope * 3600.0 * 1E3

            # later + on 23/01/2017
            y1 = (ra0-ra_fit.intercept) * 3600.0 * 1E3
            y2 = (dec0-dec_fit.intercept) * 3600.0 * 1E3

            # Red for 2MASS | later + on 23/01/2017
            ax0[row][0].scatter(x0[0], y1[0], marker='o', facecolor='r',
                                alpha=0.5, edgecolor='none')
            ax0[row][1].scatter(x0[0], y2[0], marker='o', facecolor='r',
                                alpha=0.5, edgecolor='none')

            # Blue for SDSS | later + on 23/01/2017
            ax0[row][0].scatter(x0[1:], y1[1:], marker='o', facecolor='b',
                                alpha=0.5, edgecolor='none')
            ax0[row][1].scatter(x0[1:], y2[1:], marker='o', facecolor='b',
                                alpha=0.5, edgecolor='none')

            # later + on 23/01/2017
            t_x  = np.array([-5,10])
            t_y1 = (ra_fit.slope*t_x) * 3600.0 * 1E3
            t_y2 = (dec_fit.slope*t_x) * 3600.0 * 1E3
            ax0[row][0].plot(t_x, t_y1, 'r--')
            ax0[row][1].plot(t_x, t_y2, 'r--')

            # later + on 23/01/2017
            ax0[row][0].annotate(tab2['ID'][ii], [0.05,0.95], ha='left',
                                 va='top', xycoords='axes fraction',
                                 weight='semibold')

            # later + on 23/01/2017
            s_pRA  = r'$\mu_{\alpha}$ = %+0.3f [mas/yr]' % pra0[with_2mass[ii]]
            s_pDec = r'$\mu_{\delta}$ = %+0.3f [mas/yr]' % pdec0[with_2mass[ii]]

            ax0[row][0].set_xlim(t_x)
            ax0[row][1].set_xlim(t_x)

            # later + on 23/01/2017
            ax0[row][0].annotate(s_pRA, [0.05,0.05], ha='left',
                                 va='bottom', xycoords='axes fraction')
            ax0[row][1].annotate(s_pDec, [0.05,0.05], ha='left',
                                 va='bottom', xycoords='axes fraction')

            # later + on 23/01/2017
            str_N = 'N(SDSS) = '+str(len(tab_SDSS))
            ax0[row][0].annotate(str_N, [0.95,0.95], ha='right', va='top',
                                 xycoords='axes fraction', fontsize='10')
            ax0[row][1].annotate(str_N, [0.95,0.95], ha='right', va='top',
                                 xycoords='axes fraction', fontsize='10')

            # later + on 23/01/2017
            if (ii == len2-1) and (len2 % n_panels != 0):
                for aa in np.arange(len2 % n_panels, n_panels):
                    ax0[aa][0].axis('off')
                    ax0[aa][1].axis('off')

            # later + on 23/01/2017
            if (row == n_panels-1) or (ii == len(tab2)-1):
                ax0[row][0].set_xlabel('Epoch - 2000.0')
                ax0[row][1].set_xlabel('Epoch - 2000.0')
                ax0[2][0].set_ylabel('RA - RA(J2000) [mas]') # Bug found here with index
                ax0[2][1].set_ylabel('Dec - Dec(J2000) [mas]')
                ax0[2][1].yaxis.set_label_position("right")
                fig.set_size_inches(8,8)
                fig.savefig(pp, format='pdf', bbox_inches='tight')
        #endfor
    #endif

    # + on 23/01/2017
    if silent == False: print '### Writing : ', out_pdf
    pp.close()

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
