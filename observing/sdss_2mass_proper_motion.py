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

from astropy.table import Table, Column, vstack

from astropy import coordinates as coords
from astroquery.sdss import SDSS
from astroquery.irsa import Irsa as IRSA
from astroquery.vizier import Vizier # + on 24/01/2016

from astropy.time import Time

from scipy.stats import linregress

from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust # + on 25/01/2017

SDSS_fld      = ['ra','dec','raErr','decErr','objid','run','rerun','camcol',
                 'field','obj', 'type','mode','mjd']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def linregress_slope_err(x1, fit1):
    '''
    Determine error in the slope from linear regression

    Parameters
    ----------
    x1 : np.array
      Numpy array of x values

    fit1 : scipy.stats._stats_mstats_common.LinregressResult
      Result from scipy.stats.linregress

    Returns
    -------
    sd_slope: float
      The error on the slope

    Notes
    -----
    Created by Chun Ly, 28 January 2017
     - Note that the equation below are from wikipedia:
       https://en.wikipedia.org/wiki/Regression_analysis
     - See the following page for more info:
       https://mail.scipy.org/pipermail/scipy-user/2008-May/016777.html
    '''

    mx  = x1.mean()
    sx2 = ((x1-mx)**2).sum()
    sd_slope = fit1.stderr * np.sqrt(1./sx2)
    return sd_slope
#enddef

def pm_position(tab0, epoch, silent=False, verbose=True):
    '''
    Determine RA/Dec at a given epoch

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy table from UCAC4 or MoVeRS containing J2000 and proper motion
      info. Note that ICRS WCS is adopted by default

    epoch : float
      The date in decimal year (e.g., 2017.00)

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 26 January 2017
    Modified by Chun Ly, 27 January 2017
     - Return c0 (coordinates for J2000)
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.pm_position() | '+systime()

    c0 = coords.SkyCoord(tab0['_RAJ2000'], tab0['_DEJ2000'], 'icrs',
                         unit=(u.deg))

    dtime = epoch - 2000.0 # in years

    # Note: pmRA for tab0 is in fact RA*cos(Dec). See Table 3 in:
    # http://iopscience.iop.org/article/10.1088/0004-6256/145/2/44/
    # Vizier notes say it is RA*cos(Dec) for both MoVeRS and UCAC4
    # UCAC later released results with the cos(Dec) factor
    pmRA = tab0['pmRA'] * dtime/np.cos(np.radians(c0.dec.value))
    ra_off  = coords.Angle(pmRA, unit=u.marcsec)
    dec_off = coords.Angle(tab0['pmDE'] * dtime, unit=u.marcsec)

    new_c0 = coords.SkyCoord(c0.ra + ra_off, c0.dec + dec_off, 'icrs')

    if silent == False:
        print '### End sdss_2mass_proper_motion.pm_position() | '+systime()

    return new_c0, c0
#enddef

def query_movers(c_arr, col_ID, silent=False, verbose=True):
    '''
    Query MoVeRS catalog (Theissen et al. 2016, AJ, 151, 41) to see if sample
    is present

    Parameters
    ----------
    c_arr : astropy.coordinates.sky_coordinate.SkyCoord
      Set of coordinates to cross-match against

    col_ID : astropy.table.column.Column
      Astropy column containing the ID

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 24 January 2017
    Modified by Chun Ly, 25 January 2017
     - Output full MoVeRS table instead of just an astropy.table with proper
       motion
     - Add col_ID input to include in table
    Modified by Chun Ly, 8 July 2017
     - Bug fix: Handle empty movers table
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.query_movers() | '+systime()

    n_sources = len(c_arr)

    cnt = 0
    for ii in range(n_sources):
        tab0 = Vizier.query_region(c_arr[ii], radius=5*u.arcsec,
                                   catalog='J/AJ/151/41')
        if len(tab0) != 0:
            if cnt == 0:
                # Handle case when first match is not the first source
                test_tab = tab0[0].copy()
                test_tab.add_row()
                test_tab = Table(test_tab[-1])
                for jj in range(ii):
                    test_tab.add_row()
                movers_tab = test_tab
                movers_tab = vstack([movers_tab, tab0[0]])
            else:
                movers_tab = vstack([movers_tab, tab0[0]])
            cnt += 1
        else:
            if cnt != 0: movers_tab.add_row()

    if silent == False:
        print '### End sdss_2mass_proper_motion.query_movers() | '+systime()

    # Mod on 08/07/2017
    if cnt == 0:
        ra  = np.zeros(n_sources)
        dec = np.zeros(n_sources)
        movers_tab = Table([ra,dec], names=('_RAJ2000','_DEJ2000'))
    else:
        movers_tab.add_column(col_ID, 0) # + on 25/01/2017
    return movers_tab
#enddef

def query_ucac4(c_arr, col_ID, silent=False, verbose=True):
    '''
    Query UCAC4 catalog (Zacharias et al. 2013, AJ, 145, 44) to get proper
    motion information

    Parameters
    ----------
    c_arr : astropy.coordinates.sky_coordinate.SkyCoord
      Set of coordinates to cross-match against

    col_ID : astropy.table.column.Column
      Astropy column containing the ID

    Returns
    -------
    ucac_tab : astropy.table.table.Table
      Vizier table containing all information, including proper motion

    Notes
    -----
    Created by Chun Ly, 25 January 2017
     - Add col_ID input to include in table (later added)
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.query_ucac4() | '+systime()

    n_sources = len(c_arr)

    cnt = 0
    for ii in range(n_sources):
        tab0 = Vizier.query_region(c_arr[ii], radius=5*u.arcsec,
                                   catalog='I/322A')
        if len(tab0) != 0:
            if cnt == 0:
                ucac_tab = tab0[0]
            else:
                ucac_tab = vstack([ucac_tab, tab0[0]])
            cnt += 1
        else:
            ucac_tab.add_row()
    #endfor
    if silent == False: print '## cnt : ', cnt

    if silent == False:
        print '### End sdss_2mass_proper_motion.query_ucac4() | '+systime()

    ucac_tab.add_column(col_ID, 0) # later + on 25/01/2017

    return ucac_tab
#enddef

def query_ucac5(c_arr, col_ID, silent=False, verbose=True):
    '''
    Query UCAC5 catalog (Zacharias et al. 2017, AJ, 153, 166) to get proper
    motion information

    Parameters
    ----------
    c_arr : astropy.coordinates.sky_coordinate.SkyCoord
      Set of coordinates to cross-match against

    col_ID : astropy.table.column.Column
      Astropy column containing the ID

    Returns
    -------
    ucac_tab : astropy.table.table.Table
      Vizier table containing all information, including proper motion

    Notes
    -----
    Created by Chun Ly, 9 October 2017
     - Started as copy of query_ucac4()
    Modified by Chun Ly, 23 October 2017
     - Bug fix: When first source has no UCAC catalog
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.query_ucac5() | '+systime()

    n_sources = len(c_arr)

    cnt = 0

    # + on 23/10/2017
    temp = Vizier.get_catalogs(catalog='I/340/ucac5')
    ucac_tab = Table(dtype=temp[0].dtype)

    for ii in range(n_sources):
        tab0 = Vizier.query_region(c_arr[ii], radius=5*u.arcsec,
                                   catalog='I/340/ucac5')

        if len(tab0) != 0:
            if len(ucac_tab) == 0: # Mod on 23/10/2017
                ucac_tab.add_row(tab0[0])
            else:
                ucac_tab = vstack([ucac_tab, tab0[0]])
            cnt += 1
        else:
            ucac_tab.add_row()

    #endfor
    if silent == False: print '## cnt : ', cnt

    if silent == False:
        print '### End sdss_2mass_proper_motion.query_ucac5() | '+systime()

    ucac_tab.add_column(col_ID, 0) # later + on 25/01/2017

    return ucac_tab
#enddef

def main(tab0, out_pdf=None, UCAC5=False, silent=False, verbose=True):
    '''
    Main function to determine proper motion based on 2MASS and SDSS
    coordinates

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
    Modified by Chun Ly, 24 January 2017
     - Call query_movers() and annotate panel with MoVeRS information
    Modified by Chun Ly, 25 January 2017
     - Call query_ucac4() and draw proper motion (dashed red) over 2MASS and
       SDSS data
     - Change call to query_movers() and draw proper motion (dashed green)
       over 2MASS and SDSS data
     - Output ASCII table containing proper motion
    Modified by Chun Ly, 28 January 2017
     - Overlay error bars for SDSS data
     - Get linear regression standard error on the proper motion value
     - Include cos(dec) factor for muRA
    Modified by Chun Ly, 29 January 2017
     - Pass J2000 coordinates and number of coordinate data points
    Modified by Chun Ly, 8 July 2017
     - Bug fix: Handle empty movers table
     - Overwrite table if files exist
    Modified by Chun Ly, 9 October 2017
     - Add UCAC5 keyword; handle UCAC5 case
     - Change output file suffix for UCAC5; Update annotation text for UCAC5
    '''

    if silent == False:
        print '### Begin sdss_2mass_proper_motion.main() | '+systime()

    len0 = len(tab0)

    out_path = os.path.dirname(out_pdf)+'/' # + on 25/01/2017

    with_2mass = np.array([idx for idx in range(len(tab0)) if
                           'XXXX' not in tab0['date_2mass'][idx]])

    # Initialize | + on 22/01/2017
    pra0  = np.repeat(-999.999, len0)
    pdec0 = np.repeat(-999.999, len0)

    # + on 29/01/2017
    e_pra0  = np.repeat(-999.999, len0)
    e_pdec0 = np.repeat(-999.999, len0)

    # + on 29/01/2017
    N_pm_points = np.repeat(0, len0)
    raj2000  = np.repeat(-999.99999, len0)
    decj2000 = np.repeat(-999.99999, len0)

    scale1 = 3600.0 * 1E3 # conversion to mas from deg | + on 28/01/2017

    # Deal with those without 2MASS data | + on 22/01/2017
    if len(with_2mass) > 0:
        tab2 = tab0.copy()
        tab2 = tab2[with_2mass]
        len2 = len(tab2) # later + on 23/01/2017
        c_sdss  = coords.SkyCoord(ra=tab2['ra'], dec=tab2['dec'], unit=(u.deg))
        c_2mass = coords.SkyCoord(ra=tab2['ra_2mass'], dec=tab2['dec_2mass'],
                                  unit=(u.deg))

        movers_tab = query_movers(c_sdss, tab2['ID']) # + on 24/01/2017
        print movers_tab

        # + on 25/01/2017
        if movers_tab['_RAJ2000'][0] != 0:
            outfile1 = out_path+'Proper_Motions_Alignment_Stars.MoVeRS.txt'
            if silent == False: print '### Writing : ', outfile1
            movers_tab.write(outfile1, format='ascii.fixed_width_two_line',
                             overwrite=True)

        # Mod on 09/10/2017
        if UCAC5 == False:
            ucac_tab = query_ucac4(c_sdss, tab2['ID']) # + on 25/01/2017
        else:
            ucac_tab = query_ucac5(c_sdss, tab2['ID'])
        #print ucac_tab

        # + on 25/01/2017. Mod on 09/10/2017
        outfile2 = out_path+'Proper_Motions_Alignment_Stars.UCAC4.txt'
        if UCAC5 == True:
            outfile2 = outfile2.replace('UCAC4','UCAC5')

        if silent == False: print '### Writing : ', outfile2
        ucac_tab.write(outfile2, format='ascii.fixed_width_two_line',
                       overwrite=True)

        ra_2mass  = c_2mass.ra.value
        dec_2mass = c_2mass.dec.value

        # later + on 23/01/2017
        if out_pdf == None: out_pdf = 'sdss_2mass_proper_motion.pdf'
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
            SDSS_raErr  = tab_SDSS['raErr'].data # in arcsec
            SDSS_dec    = tab_SDSS['dec'].data
            SDSS_decErr = tab_SDSS['decErr'].data # in arcsec

            t_2mass = Time(tab2['date_2mass'][ii], format='iso')
            t_2mass.format='decimalyear'

            time0 = np.array([t_2mass.value] + list(SDSS_date.value))
            ra0   = [ra_2mass[ii]]  + SDSS_ra.tolist()
            dec0  = [dec_2mass[ii]] + SDSS_dec.tolist()

            x0 = time0 - 2000.0

            # Mod on 23/01/2017
            ra_fit  = linregress(x0, ra0)
            dec_fit = linregress(x0, dec0)

            # cos(Dec) factor for muRA + on 28/01/2017
            cosd = np.cos(np.radians(c_sdss[ii].dec.value))

            # Mod on 28/01/2017
            pra0[with_2mass[ii]]  = ra_fit.slope  * scale1 * cosd # mas/yr
            pdec0[with_2mass[ii]] = dec_fit.slope * scale1        # mas/yr

            # + on 28/01/2017. Later mod to include cosd factor
            ra_fit_slope_err  = linregress_slope_err(x0, ra_fit)  * scale1*cosd
            dec_fit_slope_err = linregress_slope_err(x0, dec_fit) * scale1

            # + on 29/01/2017
            e_pra0[with_2mass[ii]]  = ra_fit_slope_err
            e_pdec0[with_2mass[ii]] = dec_fit_slope_err

            # later + on 23/01/2017
            y1 = (ra0  - ra_fit.intercept)  * scale1
            y2 = (dec0 - dec_fit.intercept) * scale1

            # + on 29/01/2017
            raj2000[with_2mass[ii]]  = ra_fit.intercept
            decj2000[with_2mass[ii]] = dec_fit.intercept

            # Red for 2MASS | later + on 23/01/2017
            ax0[row][0].scatter(x0[0], y1[0], marker='o', facecolor='r',
                                alpha=0.5, edgecolor='none')
            ax0[row][1].scatter(x0[0], y2[0], marker='o', facecolor='r',
                                alpha=0.5, edgecolor='none')

            # Blue for SDSS | later + on 23/01/2017
            # Mod on 28/01/2017 to overlay errorbars
            ax0[row][0].errorbar(x0[1:], y1[1:], yerr=SDSS_raErr*1E3,
                                 marker='o', mfc='b', alpha=0.5, capsize=0,
                                 fmt='', mec='none', linestyle='None')
            ax0[row][1].errorbar(x0[1:], y2[1:], yerr=SDSS_decErr*1E3,
                                 marker='o', mfc='b', alpha=0.5, capsize=0,
                                 fmt='', mec='none', linestyle='None')
            #ax0[row][0].scatter(x0[1:], y1[1:], marker='o', facecolor='b',
            #                    alpha=0.5, edgecolor='none')
            #ax0[row][1].scatter(x0[1:], y2[1:], marker='o', facecolor='b',
            #                    alpha=0.5, edgecolor='none')

            # later + on 23/01/2017
            t_x  = np.array([-5,10])
            t_y1 = (ra_fit.slope*t_x)  * 3600.0 * 1E3
            t_y2 = (dec_fit.slope*t_x) * 3600.0 * 1E3
            ax0[row][0].plot(t_x, t_y1, 'b--')
            ax0[row][1].plot(t_x, t_y2, 'b--')

            # Draw UCAC4 proper motion (red dashed) | + on 25/01/2017
            t_ucac = ucac_tab[ii]
            if t_ucac['_RAJ2000'] != 0.0:
                # Mod on 26/01/2017 to remove the cos(Dec) factor
                # Mod on 26/01/2017 to put the cos(Dec) factor back in
                #cosd = np.cos(np.radians(t_ucac['_DEJ2000']))
                u_pRA,  u_e_pRA  = t_ucac['pmRA'], t_ucac['e_pmRA']
                u_pDec, u_e_pDec = t_ucac['pmDE'], t_ucac['e_pmDE']

                ra_diff  = t_ucac['_RAJ2000'] - ra_fit.intercept
                dec_diff = t_ucac['_DEJ2000'] - dec_fit.intercept
                ucac_y1 = ra_diff  + u_pRA  * t_x
                ucac_y2 = dec_diff + u_pDec * t_x

                ax0[row][0].plot(t_x, ucac_y1, 'r--')
                ax0[row][1].plot(t_x, ucac_y2, 'r--')

                s_pRA =r'$\mu_{\alpha}$cos($\delta$) = %+0.3f$\pm$%0.3f' % \
                        (u_pRA,u_e_pRA)+' (UCAC4)\n'
                s_pDec=r'$\mu_{\delta}$ = %+0.3f$\pm$%0.3f' % \
                        (u_pDec,u_e_pDec)+' (UCAC4)\n'
                # + on 09/10/2017
                if UCAC5 == True:
                    s_pRA  = s_pRA.replace('UCAC4','UCAC5')
                    s_pDec = s_pDec.replace('UCAC4','UCAC5')
            else:
                s_pRA, s_pDec = '', ''

            # later + on 23/01/2017
            ax0[row][0].annotate(tab2['ID'][ii], [0.05,0.95], ha='left',
                                 va='top', xycoords='axes fraction',
                                 weight='semibold')

            # Draw MoVeRS proper motion (green dashed) | + on 24/01/2017
            # Mod on 25/01/2017
            t_movers = movers_tab[ii]
            if t_movers['_RAJ2000'] != 0.0:
                # Mod on 26/01/2017 to remove the cos(Dec) factor
                # Mod on 26/01/2017 to put the cos(Dec) factor back in
                #cosd = np.cos(np.radians(t_movers['_DEJ2000']))
                m_pRA,  m_e_pRA  = t_movers['pmRA'], t_movers['e_pmRA']
                m_pDec, m_e_pDec = t_movers['pmDE'], t_movers['e_pmDE']

                ra_diff  = t_movers['_RAJ2000'] - ra_fit.intercept
                dec_diff = t_movers['_DEJ2000'] - dec_fit.intercept
                movers_y1 = ra_diff  + m_pRA  * t_x
                movers_y2 = dec_diff + m_pDec * t_x

                ax0[row][0].plot(t_x, movers_y1, 'g--')
                ax0[row][1].plot(t_x, movers_y2, 'g--')

                s_pRA +=r'$\mu_{\alpha}$cos($\delta$) = %+0.3f$\pm$%0.3f' % \
                         (m_pRA,m_e_pRA)+' (MoVeRS)\n'
                s_pDec+=r'$\mu_{\delta}$ = %+0.3f$\pm$%0.3f' % \
                         (m_pDec,m_e_pDec)+' (MoVeRS)\n'

            # 2MASS-SDSS proper motion
            # later + on 23/01/2017
            # Mod on 28/01/2017 to report errors
            s_pRA  += r'$\mu_{\alpha}$cos($\delta$) = %+0.3f$\pm$%0.3f [mas/yr]' % \
                      (pra0[with_2mass[ii]], ra_fit_slope_err)
            s_pDec += r'$\mu_{\delta}$ = %+0.3f$\pm$%0.3f [mas/yr]' % \
                      (pdec0[with_2mass[ii]], dec_fit_slope_err)

            ax0[row][0].set_xlim(t_x)
            ax0[row][1].set_xlim(t_x)

            # later + on 23/01/2017
            ax0[row][0].annotate(s_pRA, [0.05,0.05], ha='left', fontsize=11,
                                 va='bottom', xycoords='axes fraction')
            ax0[row][1].annotate(s_pDec, [0.05,0.05], ha='left', fontsize=11,
                                 va='bottom', xycoords='axes fraction')

            N_pm_points[with_2mass[ii]] = len(tab_SDSS)+1 # + on 29/01/2017

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
                subplots_adjust(left=0.10, right=0.95, top=0.98, bottom=0.10,
                                hspace=0.09) # + on 25/01/2017
                fig.set_size_inches(8,8)
                fig.savefig(pp, format='pdf') #, bbox_inches='tight')
            else: # Mod on 25/01/2017 to avoid x-axis tick labels
                ax0[row][0].set_xticklabels([])
                ax0[row][1].set_xticklabels([])
        #endfor
    #endif

    # + on 23/01/2017
    if silent == False: print '### Writing : ', out_pdf
    pp.close()

    if silent == False:
        print '### End sdss_2mass_proper_motion.main() | '+systime()

    # Mod on 29/01/2017
    names0 = ('ID', '_RAJ2000', '_DEJ2000', 'pmRA', 'e_pmRA', 'pmDE',
              'e_pmDE', 'N_pm_points')
    vec0   = [tab0['ID'], raj2000, decj2000, pra0, e_pra0,
              pdec0, e_pdec0, N_pm_points]
    sdss_2mass_tab = Table(vec0, names=names0)
    return sdss_2mass_tab, movers_tab, ucac_tab
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
