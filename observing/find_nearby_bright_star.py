"""
find_nearby_bright_star
=======================

Search SDSS or 2MASS catalogs for nearest bright star as offset star
or to determine PA for longslit spectroscopy. Also generates
finding charts

Requirements:
 astroquery (https://astroquery.readthedocs.io/en/latest/)
 aplpy (https://aplpy.github.io/)
 pdfmerge ('pip install pdfmerge' will install it)
"""

import sys, os

from chun_codes import systime, match_nosort_str
import string

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
from astroquery.skyview import SkyView
from astropy.time import Time # + on 21/01/2017

import sdss_2mass_proper_motion as pm # in Zcalbase_gal.observing
import aplpy # + on 24/12/2016

from PyMontage.scripts import montage_reproj # + on 02/01/2017

# For SDSS only
SDSS_fld      = ['ra','dec','objid','run','rerun','camcol','field','obj',
                 'type','mode','mjd']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def get_PA(c0, c1, slitlength=99*u.arcsec, MMT=False, silent=True, verbose=False):
    '''
    Function to determine the PA on the sky and central coordinates given
    two WCS coordinates

    Parameters
    ----------
    c0 : `astropy.coordinates` object
      Central coordinate of target

    c1 : `astropy.coordinates` object
      Central coordinate of nearby bright star

    MMT : boolean
      Indicate if finding charts are for MMT (makes it simpler)

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------
    PA : float
      Position angle between c0 and c1. Positive values is for East of North

    c_ctr : `astropy.coordinates` object
      Coordinate associated with center between target and nearby bright star

    Notes
    -----
    Created by Chun Ly, 3 January 2017
     - Added slitlength option to get coordinates of longslit
    Modified by Chun Ly, 4 January 2017
     - Get PA for longslit for illustration on finding chart. This PA is
       different from SkyCoords.position_angle(). This is merely to show the
       slit length
    Modified by Chun Ly, 9 January 2017
     - Fix longslit coordinates. Correct ra0 and dec0 this time
    Modified by Chun Ly, 24 January 2017
     - Force fk5 coordinate for astropy.coords
    Modified by Chun Ly,  8 October 2017
     - Handle MMT: Use bright star coordinate as center
    '''
    
    if silent == False:
        print '### Begin find_nearby_bright_star.get_PA | '+systime()

    PA  = c0.position_angle(c1).degree # + => East of North
    PAr = np.radians(PA)

    # Get central coordinate
    # Mod on 08/10/2017
    if MMT == False:
        ra_avg  = np.average([c0.ra.value, c1.ra.value])
        dec_avg = np.average([c0.dec.value, c1.dec.value])
    else:
        ra_avg, dec_avg = c1.ra.value, c1.dec.value

    # Mod on 24/01/2017 for fk5
    c_ctr = coords.SkyCoord(ra_avg, dec_avg, 'fk5', unit=(u.deg))

    # Get edges of longslit | Added later
    # Mod on 04/01/2017
    # t_y, t_x = c_ctr.dec.degree - c0.dec.degree, c_ctr.ra.degree - c0.ra.degree
    # This is different from PA. 90 deg is because of arctan2 definition.
    # PA2=0 correspond to N
    # PA2  = 90 - np.degrees(np.arctan2(t_y, t_x)) # in degree.
    # PA2r = np.radians(PA2)
    ra0  = 0.5*slitlength.to(u.arcsec).value * np.sin(PAr) / np.cos(np.radians(dec_avg))
    dec0 = 0.5*slitlength.to(u.arcsec).value * np.cos(PAr)

    ra_offset  = coords.Angle(ra0, unit=u.arcsec)
    dec_offset = coords.Angle(dec0, unit=u.arcsec)

    new_pos = coords.SkyCoord(c_ctr.ra+ra_offset*[-1,1],
                              c_ctr.dec+dec_offset*[-1,1])

    longslit_list = []
    longslit_list.append(np.array([[new_pos[0].ra.value, new_pos[1].ra.value],
                                   [new_pos[0].dec.value, new_pos[1].dec.value]]))

    if silent == False:
        print '### End find_nearby_bright_star.get_PA | '+systime()

    return PA, c_ctr, longslit_list
#enddef

def get_offsets(c_ref, c0):
    '''
    Get relative offsets between two positions and provide values in units
    of arcsec

    Parameters
    ----------
    c_ref : `astropy.coordinates` object
      Central coordinate of reference (i.e., the position to offset FROM)

    c0 : `astropy.coordinates` object
      Central coordinate of target (i.e., the position to offset TO)

    Returns
    -------
    dra0 : float
      Right ascension offsets in arcsec

    ddec0 : float
      Declination offsets in arcsec

    dist0 : float
      Distance between c0 and c_ref in arcsec

    Notes
    -----
    Created by Chun Ly, 10 January 2017
    Modified by Chun Ly, 27 January 2017
     - Check to make sure frames are equivalent
    Modified by Chun Ly, 15 February 2017
     - Bug found with spherical_offsets_to() for large offsets.
       Re-computing using standard definition:
        d(RA)  = (RA - RA0)  * cos(DEC0)
        d(Dec) = (Dec - Dec0)
    '''

    # + on 27/01/2017 for frame check
    if c_ref.is_equivalent_frame(c0) == False:
        c_ref0 = c_ref.transform_to(c0.frame)
    else: c_ref0 = c_ref

    # Give coordinate offsets from reference [c_ref0] to target [c0]
    # Mod on 15/02/2017
    dra0  = (c0.ra.deg - c_ref0.ra.deg)   * 3600.0 * np.cos(c_ref0.dec.radian)
    ddec0 = (c0.dec.deg - c_ref0.dec.deg) * 3600.0
    # dra0, ddec0 = c_ref0.spherical_offsets_to(c0)
    # dra0  = dra0.to(u.arcsec).value
    # ddec0 = ddec0.to(u.arcsec).value
    dist0 = np.sqrt(dra0**2+ddec0**2)
    return dra0, ddec0, dist0
#enddef

def plot_finding_chart(fitsfile, t_ID, band0, c0, c1, mag_str, mag_table,
                       out_pdf=None, slitlength=99*u.arcsec, catalog='SDSS',
                       image=None, pmfix=False, c1_new=None, c1_2000=None,
                       epoch=2000.0, MMT=False, silent=True, verbose=False):
    '''
    Function to plot FITS images with WCS on the x- and y-axes

    Parameters
    ----------
    fitsfile : string
      Filename of FITS file

    t_ID : string
      Name of source

    band0 : string
      Observation waveband of image

    c0 : `astropy.coordinates` object
      Central coordinate of target

    c1 : `astropy.coordinates` object
      Central coordinate of nearby bright stars

    mag_str : string
      Information on SDSS ugriz and 2MASS JHK for target

    mag_table : astropy.table.table.Table
      Table containing SDSS and 2MASS information, and proper motion

    out_pdf : string
      Output PDF filename.
      Default: based on [fitsfile], replacing '.fits.gz' or '.fits' with '.pdf'

    MMT : boolean
      Indicate if finding charts are for MMT (makes it simpler)

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    catalog : string
      Which survey (e.g., SDSS, 2MASS) the finding chart image is from.
      Default: 'SDSS'
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 24 December 2016
    Modified by Chun Ly, 2 January 2017
     - Fix previous bug with show_markers (change facecolor)
     - Added input of c1 to indicate where bright alignment stars are
     - Added t_ID, band0 inputs for plt.annotate call
     - Changed out_pdf to keyword rather than input variable
     - Added catalog keyword
     - Change axis label values to clean it up
    Modified by Chun Ly, 3 January 2017
     - Call get_PA() and draw slit
     - Add annotated text in the lower left corner for info about
       long-slit alignment
    Modified by Chun Ly, 4 January 2017
     - Added slitlength keyword input to pass to get_PA()
    Modified by Chun Ly, 7 January 2017
     - Change major tick labels for Dec
    Modified by Chun Ly, 8 January 2017
     - Add offset values to central position
     - Add c1_mag input to get magnitudes for bright adjacent stars
    Modified by Chun Ly, 9 January 2017
     - Switch [c1_mag] to string of mag from sdss_mag_str(), [mag_str]
    Modified by Chun Ly, 10 January 2017
     - Call get_offsets() to simply get values
    Modified by Chun Ly, 19 January 2017
     - Add image keyword option
    Modified by Chun Ly, 22 January 2017
     - Add [mag_table] array input to get proper motion
     - Overlay 2MASS positions in red
    Modified by Chun Ly, 23 January 2017
     - Add SDSS observation date to label
    Modified by Chun Ly, 26 January 2017
     - Add pmfix keyword option
    Modified by Chun Ly, 27 January 2017
     - Add c1_new, c1_2000, and epoch keywords option to update coordinates
       and annotated text
    Modified by Chun Ly, 28 January 2017
     - Add source of proper motion to [bt_txt]
     - Improve annotated text for pmfix case and when image overlay is 2MASS
    Modified by Chun Ly, 29 January 2017
     - Change annotated text for proper motion if c1_new is not available
    Modified by Chun Ly, 8 October 2017
     - Add MMT keyword
     - Pass MMT keyword to get_PA()
    Modified by Chun Ly, 9 October 2017
     - Avoid plotting blue circles for pmfix == True
     - Adjust text label for MMT
    Modified by Chun Ly, 10 October 2017
     - Handle obs date for 2MASS catalog
    '''

    # + on 19/01/2017
    if image == None:
        if catalog == 'SDSS':
            band_str0 = mag_filt.replace('modelMag_','')
        if catalog == '2MASS':
            band_str0 = mag_filt.replace('_m','').upper()
        image = catalog+'-'+band_str0
        if silent == False: '## image : ', image

    # + on 02/01/2017
    if out_pdf == None:
        out_pdf = fitsfile.replace('.fits.gz', '.pdf').replace('.fits','.pdf')

    gc = aplpy.FITSFigure(fitsfile, figsize=(8,8), north=True)

    gc.show_grayscale(invert=True)
    gc.set_tick_color('black')

    # + on 02/01/2017 to fix plotting
    # Mod on 07/01/2017
    gc.set_tick_yspacing('auto')
    gc.ticks.set_yspacing(1/60.0) # Every 1 arcmin in Dec
    gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')
    
    # Draw slit between target and nearest bright star | + on 03/01/2017
    # Mod on 27/01/2017 to handle pmfix case
    do_pm = False
    if pmfix == True and c1_new != None:
        do_pm = True
        print '## Using proper-motion based position'
        PA, c_ctr, longslit_list = get_PA(c0, c1_new, slitlength=slitlength,
                                          MMT=MMT)
    else:
        PA, c_ctr, longslit_list = get_PA(c0, c1[0], slitlength=slitlength,
                                          MMT=MMT)

    gc.show_lines(longslit_list, layer='slit', color='black', linewidth=0.5)
    #gc.show_rectangles([c_ctr.ra.value], [c_ctr.dec.value], 1/3600.0, 99/3600.0)

    # Fix bug. marker='+' won't work with facecolor='none'
    gc.show_markers([c0.ra.value], [c0.dec.value], layer='primary',
                    edgecolor='red', facecolor='red', marker='+', s=25)

    # + on 02/01/2017
    # Mod on 09/10/2017
    if pmfix != True:
        gc.show_markers([c1.ra.value], [c1.dec.value], layer='secondary',
                        edgecolor='blue', facecolor='none', marker='o', s=25,
                        linewidth=0.5)

    # Draw green circle for new position at epoch | + on 27/01/2017
    if do_pm == True:
        gc.show_markers([c1_new.ra.value], [c1_new.dec.value], layer='',
                        edgecolor='green', facecolor='none', marker='o', s=25,
                        linewidth=0.5)

    # Overlay 2MASS coordinates in red | + on 22/01/2017
    if not MMT: # Mod on 08/10/2017
        gc.show_markers(mag_table['ra_2mass'], mag_table['dec_2mass'],
                        edgecolor='red', facecolor='none', marker='o',
                        s=15, linewidth=0.5)
    
    # Label things in lower left text | + on 03/01/2017
    str_c_t  = c0.to_string('hmsdms').split(' ')
    str_c    = c_ctr.to_string('hmsdms').split(' ')
    str_c_bt = c1[0].to_string('hmsdms').split(' ')

    # Mod on 27/01/2017 for pmfix case
    if do_pm == False:
        dra, ddec, dist    = get_offsets(c1[0], c0)    # Mod on 10/01/2017
        dra2, ddec2, dist2 = get_offsets(c1[0], c_ctr) # Mod on 10/01/2017
    else:
        dra, ddec, dist    = get_offsets(c1_new, c0)    # Mod on 10/01/2017
        dra2, ddec2, dist2 = get_offsets(c1_new, c_ctr) # Mod on 10/01/2017

    # + on 27/01/2017
    if do_pm == True:
        str_c_bt_new  = c1_new.to_string('hmsdms').split(' ')
        str_c_bt_2000 = c1_2000.to_string('hmsdms').split(' ')

    # + on 23/01/2017. Mod on 10/10/2017
    if catalog == 'SDSS':
        obs_date = Time(mag_table['mjd'][0], format='mjd')
    if catalog == '2MASS':
        obs_date = Time(mag_table['xdate'][0])
    obs_date.format='iso'

    # Mod on 09/01/2017, 10/01/2017, 22/01/2017, 27/01/2017, 28/01/2017
    bt_txt  = 'Target: RA='+str_c_t[0]+', Dec='+str_c_t[1]+'\n\n'

    if not MMT:
        bt_txt += 'Slit Center: RA='+str_c[0]+', Dec='+str_c[1]
    else:
        bt_txt += 'Slit Align. Coord.: RA='+str_c[0]+', Dec='+str_c[1]

    if do_pm == True:
        bt_txt += '  Epoch='+str(epoch)+'\n'
    else: bt_txt += '\n'
    bt_txt += ('Slit PA = %7.3f' % PA) + ' deg'
    if do_pm == True:
        bt_txt += '  Epoch='+str(epoch)+'\n'
    else: bt_txt += '\n'

    if not MMT:
        bt_txt += 'Offsets : (%+.3f", %+.3f");  %.2f"\n\n' % (dra2, ddec2, dist2)
    else:
        bt_txt += '\n'

    if do_pm == False:
        bt_txt += 'Offset Star: RA='+str_c_bt[0]+', Dec='+str_c_bt[1]
        bt_txt += ' Date='+obs_date.value.split(" ")[0]+'\n' # + on 23/01/2017
    else:
        bt_txt += 'Offset Star: RA='+str_c_bt_2000[0]+', Dec='+str_c_bt_2000[1]
        bt_txt += '  Epoch=J2000.0\n'
        bt_txt += 'Offset Star: RA='+str_c_bt_new[0]+', Dec='+str_c_bt_new[1]
        bt_txt += '  Epoch='+str(epoch)+'\n'

    if not MMT:
        bt_txt += 'Offsets : (%+.3f", %+.3f");  %.2f"\n' % (dra, ddec, dist)
    else:
        bt_txt += 'Offset Dist to Target: (%+.3f", %+.3f");  %.2f"\n' % (dra, ddec, dist)

    if do_pm == True: # Mod on 30/01/2017
        bt_txt += r'$\mu(\alpha)$cos($\delta$) = %+.2f mas/yr, ' % mag_table['pRA'][0]
        bt_txt += r'$\mu(\delta)$ = %+.2f mas/yr, Source: %s' % \
                  (mag_table['pDec'][0], mag_table['p_source'][0])
    bt_txt += '\n'+mag_str[0]

    ctype = 'white' if ('2MASS' in fitsfile) else 'magenta' # + on 28/01/2017
    gc.add_label(0.03, 0.15, bt_txt, color=ctype, relative=True,
                 ha='left', va='bottom', weight='medium', size='small')

    # Label upper left source, catalog, and band | + on 02/01/2017
    # Mod on 19/01/2017
    lab0 = t_ID+'\nCatalog: '+catalog+'\nImage: '+image.replace('-',' ')
    gc.add_label(0.03, 0.925, lab0, relative=True, ha='left', va='top',
                 weight='bold', size='large')

    gc.savefig(out_pdf)
#enddef

def get_sdss_images(c0, out_fits, band=u'i', silent=True, verbose=False):
    '''
    Function to grab SDSS FITS image that is associated with the provided
    coordinate

    Parameters
    ----------
    c0 : str or `astropy.coordinates` object
      The target around which to search. It may be specified as a string
      in which case it is resolved using online services or as the
      appropriate `astropy.coordinates` object. ICRS coordinates may also
      be entered as strings as specified in the `astropy.coordinates`
      module.

    band : string
      Filter name for image to obtain from query (e.g., 'u', 'g', 'r', 'i', 'z')
      Default: 'i'
    
    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False
	  
    Returns
    -------
    t_hdu : astropy.io.fits.image.PrimaryHDU
      An HDU object containing FITS header and FITS data for the 2MASS
      finding chart image

    Notes
    -----
    Created by Chun Ly, 24 December 2016
    '''

    imgs = SDSS.get_images(coordinates=c0, radius=1*u.arcmin, band=band,
                           timeout=180)

    n_frames = len(imgs)

    if silent == False: print '### Writing : ', out_fits
    for ff in range(n_frames):
        t_hdu = imgs[ff][0]
        if ff == 0:
            t_hdu.writeto(out_fits, clobber=True)
        else:
            fits.append(out_fits, t_hdu.data, t_hdu.header)
    return t_hdu
#enddef

def get_2mass_images(c0, out_fits, band=u'H', silent=True, verbose=False):
    '''
    Function to grab 2MASS FITS image that is associated with the provided
    coordinate

    Parameters
    ----------
    c0 : str or `astropy.coordinates` object
      The target around which to search. It may be specified as a string
      in which case it is resolved using online services or as the
      appropriate `astropy.coordinates` object. ICRS coordinates may also
      be entered as strings as specified in the `astropy.coordinates`
      module.

    band : string
      Filter name for image to obtain from query (e.g., 'J', 'H', 'K')
      Default: 'H'

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------
    t_hdu : astropy.io.fits.image.PrimaryHDU
      An HDU object containing FITS header and FITS data for the 2MASS
      finding chart image

    Notes
    -----
    Created by Chun Ly, 19 January 2017
    Modified by Chun Ly, 20 January 2017
     - Fix bug with SkyView.get_images() call
    '''

    if silent == False:
        print '### Begin find_nearby_bright_star.get_2mass_images '+systime()

    imgs = SkyView.get_images(position=c0, survey=['2MASS-'+band],
                              width=5*u.arcmin, height=5*u.arcmin)

    if silent == False: print '### Writing : ', out_fits

    # Mod on 20/01/2017 to clean it up
    t_hdu = imgs[0][0]
    t_hdu.writeto(out_fits, clobber=True)

    if silent == False:
        print '### End find_nearby_bright_star.get_2mass_images '+systime()

    return t_hdu
#enddef

def sdss_mag_str(table, TWOMASS=True, runall=True):
    '''
    Function to create a string of magnitudes for labeling purposes on
    finding charts

    Parameters
    ----------
    table : astropy.table
      SDSS table from SDSS.query_region()

    TWOMASS : bool
      Set to True to perform cross-matching to get 2MASS photometry.
      Default: True

    Returns
    -------
    mag_str : string

    Notes
    -----
    Created by Chun Ly, 9 January 2017
    Modified by Chun Ly, 10 January 2017
     - Add option to get 2MASS photometry by running IRSA.query_region()
     - Combine SDSS and 2MASS photometry, output astropy.table
    Modified by Chun Ly, 21 January 2017
     - Added columns of 2MASS's RA, Dec, and observation date
    Modified by Chun Ly, 22 January 2017
     - Improve efficiency for runall=True with larger radius in IRSA.query_region
    Modified by Chun Ly, 24 January 2017
     - Force fk5 coordinate for astropy.coords
    Modified by Chun Ly, 28 January 2017
     - Add p_source to indicate where proper motion comes from. This will be
       updated later for UCAC4 or MoVeRS proper motion (see main())
    Modified by Chun Ly, 29 January 2017
     - Add proper motion error columns to astropy Table
    '''

    n_sources = len(table)

    # + on 10/01/2017
    mag_table = table.copy()
    mag_table.remove_columns(['run','rerun','camcol','field','obj'])

    cols0  = [a for a in table.colnames if 'modelMag_' in a]
    f_arr0 = []
    for cc in range(len(cols0)):
        t_filt = cols0[cc].replace('modelMag_','')
        cmd = "mag_"+t_filt+\
              " = ['"+t_filt+"='+('%5.2f  ' % a) for a in table[cols0["+str(cc)+"]]]"
        exec(cmd)
        f_arr0.append("mag_"+t_filt)

    q_str = string.lowercase[:len(f_arr0)]
    cmd0  = "mag_str = ["+'+'.join(q_str)+' for '+\
            ','.join(q_str)+' in zip('+','.join(f_arr0)+')]'
    exec(cmd0)

    # Get 2MASS photometry using a astroquery approach
    # + on 10/01/2017
    col_Jmag  = Column(np.repeat(-99.00, n_sources), name='j_m')
    col_Hmag  = Column(np.repeat(-99.00, n_sources), name='h_m')
    col_Kmag  = Column(np.repeat(-99.00, n_sources), name='k_m')
    col_eJmag = Column(np.repeat(-99.00, n_sources), name='j_cmsig')
    col_eHmag = Column(np.repeat(-99.00, n_sources), name='h_cmsig')
    col_eKmag = Column(np.repeat(-99.00, n_sources), name='k_cmsig')

    # + on 21/01/2017
    col_ra   = Column(np.zeros(n_sources), name='ra_2mass')
    col_dec  = Column(np.zeros(n_sources), name='dec_2mass')
    col_date = Column(np.repeat('XXXX-XX-XX', n_sources), name='date_2mass')

    if TWOMASS == True:
        # Mod on 24/01/2017 for fk5
        t_c0 = coords.SkyCoord(table['ra'], table['dec'], 'fk5', unit=(u.deg))

        # Mod on 22/01/2017 for efficiency with runall == True
        rad0 = 5*u.arcmin if runall == True else 1*u.arcsec
        if n_sources == 1: rad0 = 1*u.arcsec
        table_2mass = IRSA.query_region(t_c0[0], catalog='fp_psc',
                                        radius=rad0)

        JHK_str0 = ['JHK=None'] * n_sources
        if len(table_2mass) > 0:
            # Mod on 24/01/2017 for fk5
            # Mod on 06/10/2017 - Use full coordinate instead of rounded values
            c_2mass = coords.SkyCoord(table_2mass['clon'], table_2mass['clat'])
            #c_2mass = coords.SkyCoord(table_2mass['ra'], table_2mass['dec'],
            #                          'fk5', unit=(u.deg))

            idx_arr, idx_ref, d2d, d3d = t_c0.search_around_sky(c_2mass, 1*u.arcsec)

            for jj in range(len(idx_arr)):
                i1 = idx_arr[jj]
                i2 = idx_ref[jj]
                tab0 = table_2mass[i1]
                JHK_str0[i2] = 'J=%5.2f  H=%5.2f  K=%5.2f' % \
                               (tab0['j_m'],tab0['h_m'],tab0['k_m'])

                # + on 10/01/2017
                col_Jmag[i2]  = tab0['j_m']
                col_Hmag[i2]  = tab0['h_m']
                col_Kmag[i2]  = tab0['k_m']
                col_eJmag[i2] = tab0['j_cmsig']
                col_eHmag[i2] = tab0['h_cmsig']
                col_eKmag[i2] = tab0['k_cmsig']

                # + on 21/01/2017
                col_ra[i2]    = tab0['ra']
                col_dec[i2]   = tab0['dec']
                col_date[i2]  = tab0['xdate']
            #endfor
        #endif

        # + on 22/01/2017
        mag_str = [a + b for a,b in zip(mag_str,JHK_str0)]

        # + on 10/01/2017
        n_cols = len(mag_table.colnames)
        # Mod on 22/01/2017
        vec0 = [col_ra, col_dec, col_Jmag, col_eJmag, col_Hmag,
                col_eHmag, col_Kmag, col_eKmag, col_date]
        mag_table.add_columns(vec0, np.repeat(n_cols-1, len(vec0)))

        # + on 22/01/2017
        pRA  = np.zeros(len(mag_table))
        pDec = np.zeros(len(mag_table))

        # + on 29/01/2017
        e_pRA  = np.zeros(len(mag_table))
        e_pDec = np.zeros(len(mag_table))

        # Mod on 22/01/2017
        if len(table_2mass) > 0:
            if len(idx_arr) > 0:
                pRA, pDec = pm.old(mag_table)

        col_pRA  = Column(pRA,  name='pRA')
        col_pDec = Column(pDec, name='pDec')

        # + on 29/01/2017
        col_e_pRA  = Column(e_pRA,  name='e_pRA')
        col_e_pDec = Column(e_pDec, name='e_pDec')

        # + on 28/01/2017
        col_p_src = Column(np.repeat('SDSS-2MASS', len(pRA)), name='p_source')

        n_cols = len(mag_table.colnames)
        vec0 = [col_pRA, col_e_pRA, col_pDec, col_e_pDec, col_p_src]
        mag_table.add_columns(vec0, np.repeat(n_cols-1, len(vec0)))
    #endif

    return mag_str, mag_table
#enddef

def main(infile, out_path, finding_chart_path, finding_chart_fits_path,
         max_radius=60*u.arcsec, mag_limit=20.0, mag_filt='modelMag_i',
         catalog='SDSS', image=None, format0='commented_header',
         slitlength=99*u.arcsec, runall=True, alignment_file='', pmfix=False,
         epoch=2000.0, pm_out_pdf=None, sig_min=3.0, MMT=False, UCAC5=False,
         silent=False, verbose=True):

    '''
    Main function to find nearby star

    Parameters
    ----------
    infile : string
      Filename for ASCII file that contains ID, RA, and DEC
      RA and DEC must be provided in sexigesimal format.
      The columns should be referred to as 'ID', 'RA', and 'DEC'

    out_path : string
      Full path to output ASCII tables. Must end with a '/'

    finding_chart_path : string
      Full path for outputted PDF finding charts. Must end with a '/'

    finding_chart_fits_path : string
      Full path for outputted finding charts. Must end with a '/'

    max_radius : float
      Maximum radius. Provide with astropy.units for arcsec,
      arcmin, or degrees. Default: 60 * u.arcsec

    mag_limit : float
      Faintest source to consider in AB mag. Default: 20.0 mag

    mag_filt : string
      The name of the filter adopt the above [mag_limit]
      Default: 'modelMag_i'

    catalog : string
      The survey to extract a catalog.
      Either 'SDSS' or '2MASS'. Default: 'SDSS'
    
    image : string
      The survey to extract the finding chart image.
      Either 'SDSS' or '2MASS'. Default: None

    format : string
      Format of infile ASCII file. Default: "commented_header"

    sig_min : float
      Minimum number of sigma to use proper motion from UCAC4 or MoVERS
      Default: 3.0

    MMT : boolean
      Indicate if finding charts are for MMT (makes it simpler)

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 23 December 2016
    Modified by Chun Ly, 24 December 2016
     - Added query for 2MASS
     - Keep those with mode = 1 (Primary sources only, reduces duplication)
    Modified by Chun Ly, 02 January 2017
     - if statement if [good] is empty
     - Added different path for FITS finding chart [finding_chart_fits_path]
    Modified by Chun Ly, 03 January 2017
     - Re-define c1 so that order is consistent with [xid] sorting
    Modified by Chun Ly, 4 January 2017
     - Added slitlength keyword input to pass to plot_finding_chart()
    Modified by Chun Ly, 6 January 2017
     - Call montage_reproj() to perform image re-projection of SDSS images,
       generate new mosaic with default image sizes
     - Change input to plot_finding_chart() for new mosaicked image
    Modified by Chun Ly, 7 January 2017
     - Added check to avoid re-creating mosaic FITS image for finding chart
    Modified by Chun Ly, 8 January 2017
     - Determine [c1_mag], nearby star brightness
    Modified by Chun Ly, 9 January 2017
     - Call sdss_mag_str()
     - Re-organize to handle SDSS and 2MASS cases
    Modified by Chun Ly, 10 January 2017
     - Get table of magnitudes and write to ASCII file in out_path
    Modified by Chun Ly, 19 January 2017
     - Add call to get_2mass_images()
     - Add image keyword option
    Modified by Chun Ly, 20 January 2017
     - Handle case of SDSS catalog but with 2MASS images for finding charts
    Modified by Chun Ly, 24 January 2017
     - Force fk5 coordinate for astropy.coords
    Modified by Chun Ly, 25 January 2017
     - Add pmfix and alignment_file keyword option to read in proper motion
       information
    Modified by Chun Ly, 26 January 2017
     - Re-organize so for loop over sources works for pmfix=True
     - Call sdss_2mass_proper_motion.main()
     - Add epoch keyword to pass to sdss_2mass_proper_motion.pm_position()
     - Add pm_out_pdf keyword to pass to sdss_2mass_proper_motion.main()
    Modified by Chun Ly, 27 January 2017
     - Get J2000 coordinates from sdss_2mass_proper_motion.pm_position()
    Modified by Chun Ly, 28 January 2017
     - Update proper motion to include UCAC or MoVeRS
    Modified by Chun Ly, 29 January 2017
     - Handle 2MASS-SDSS proper motion and pass those if reliable based on
       4-sigma criteria. Adopt MoVeRS and UCAC4 first
    Modified by Chun Ly, 29 January 2017
     - Change limit to 3-sigma to include DEEP14 proper motion
    Modified by Chun Ly, 7 July 2017
     - Fix bug with ascii format
     - Fix bug when movers table is empty
    Modified by Chun Ly, 8 July 2017
     - Add sig_min keyword to allow user to specify minimum for using PM info
    Modified by Chun Ly, 8 October 2017
     - Add MMT keyword
    Modified by Chun Ly, 9 October 2017
     - Add UCAC5 keyword; Handle UCAC5 case
    Modified by Chun Ly, 10 October 2017
     - Handle using 2MASS catalog
     - Add PM info to table for both SDSS and 2MASS catalogs
    '''

    if silent == False:
        print '### Begin find_nearby_bright_star.main | '+systime()

    if catalog == 'SDSS':  band0 = mag_filt.replace('modelMag_','')
    if catalog == '2MASS': band0 = mag_filt.replace('_m','').upper()

    # + on 19/01/2017
    if image == None: image = catalog + '-' + band0
    if silent == False: print '### image : ', image

    if silent == False: print '### Reading : ', infile
    data0 = asc.read(infile, format=format0)
    ID    = data0['ID'].data
    ID0   = [ss.replace('*','') for ss in ID]

    # Cross-match against those preferred sources to get proper motion info
    # + on 26/01/2017
    if pmfix == True:
        if silent == False: print '### Reading : ', alignment_file
        a_tab0 = asc.read(alignment_file)
        a_ID0  = [str0.replace('_off','') for str0 in a_tab0['ID']]
        idx1, idx2 = match_nosort_str(ID0, a_ID0)
        print '## idx1 : ', len(idx1)
        data0 = data0[idx1]

        # SDSS-2MASS, MoVeRS, and UCAC4 proper motion catalogs
        # Mod on 09/10/2017
        t_s2, t_movers, t_ucac = pm.main(a_tab0, pm_out_pdf, UCAC5=UCAC5)

        # Adopt a 4-sigma criteria for trusting proper motion
        # Mod on 30/01/2017 to adopt 3-sigma instead
        # Mod on 07/07/2017
        print '### Using proper motion that is more reliable than '+str(sig_min)+'sigma'

        if t_movers['_RAJ2000'][0] != 0:
            m_idx = np.where((abs(t_movers['pmRA']/t_movers['e_pmRA']) >= sig_min) &
                             (abs(t_movers['pmDE']/t_movers['e_pmDE']) >= sig_min))[0]
        else: m_idx = []
        u_idx = np.where((abs(t_ucac['pmRA']/t_ucac['e_pmRA']) >= sig_min) &
                         (abs(t_ucac['pmDE']/t_ucac['e_pmDE']) >= sig_min))[0]

        # + on 29/01/2017
        pm_source = np.zeros(len(data0))
        pm_source[m_idx] = 1
        pm_source[u_idx] = 2

        # 2MASS-SDSS | + on 29/01/2017
        s2_idx = np.where((abs(t_s2['pmRA']/t_s2['e_pmRA']) >= 4.0) &
                          (abs(t_s2['pmDE']/t_s2['e_pmDE']) >= 4.0) &
                          (pm_source == 0) & (t_s2['N_pm_points'] >= 5))[0]

        if len(m_idx) > 0:
            print '### Will use MoVeRS proper motion for the following : '
            print '### : ', t_movers['ID'][m_idx]
        if len(u_idx) > 0:
            print '### Will use UCAC4 proper motion for the following : '
            print '### : ', t_ucac['ID'][u_idx]
        if len(s2_idx) > 0: # + on 29/01/2017
            print '### Will use SDSS-2MASS proper motion for the following : '
            print '### : ', t_s2['ID'][s2_idx]

        if silent == False: print '### Adopted epoch : ', epoch
        # Mod on 27/01/2017
        if len(m_idx) > 0: # Mod on 07/07/2017
            m_c0, m_c0_2000   = pm.pm_position(t_movers, epoch)
        u_c0, u_c0_2000   = pm.pm_position(t_ucac,   epoch)
        s2_c0, s2_c0_2000 = pm.pm_position(t_s2,     epoch) # + on 29/01/2017

    ID  = data0['ID'].data
    RA  = data0['RA'].data
    DEC = data0['DEC'].data
    n_sources = len(ID)

    ID0 = [ss.replace('*','') for ss in ID] # Clean things up 24/12/2016

    if pmfix == False and silent == False:
        print '## Search criteria : '
        print '## max_radius(arcsec) : ', max_radius.to(u.arcsec).value
        print '## mag_limit : ', mag_limit
        print '## filter selection : ', mag_filt

    # Mod on 10/10/2017
    mag_table0 = [] #None # Initialize | + on 10/01/2017
    for ii in range(n_sources):
        # Mod on 24/01/2017 for fk5
        c0 = coords.SkyCoord(RA[ii], DEC[ii], 'fk5', unit=(u.hour, u.degree))
        # Moved up on 09/01/2017
        out_table_file = out_path + ID0[ii] + '.'+catalog+'.nearby.txt'

        skip0 = 0

        if pmfix == True:
            if silent == False: print '### Reading : ', out_table_file
            xid = asc.read(out_table_file)
        else:
            if catalog == 'SDSS':
                xid = SDSS.query_region(c0, radius=max_radius, data_release=12,
                                        photoobj_fields=SDSS_phot_fld)

            # + on 24/12/2016 | Moved up on 09/01/2017
            if catalog == '2MASS':
                xid = IRSA.query_region(c0, catalog='fp_psc', radius=max_radius)
                good = np.where(xid[mag_filt] <= mag_limit)[0]

            # Get distance from target
            if len(xid) > 0:
                # Mod on 24/01/2017 for fk5
                c1 = coords.SkyCoord(xid['ra'], xid['dec'], 'fk5', unit=(u.deg))
                sep = c0.separation(c1).to(u.arcsec).value
                col1 = Column(sep, name='Dist(arcsec)')
                xid.add_column(col1)

            # Keep stars only (type == 6)
            # http://www.sdss.org/dr12/algorithms/classify/#photo_class
            # Keep primary target to deal with duplicate entries | 24/12/2016
            # Avoid unusual mag values | Mod on 24/12/2016
            # Fix SDSS.query_region problem with max_radius | Mod on 08/01/2017
            if catalog == 'SDSS':
                good = np.where((xid[mag_filt] <= mag_limit) &
                                (xid[mag_filt] != -9999.0) & # Mod on 24/12/2016
                                (xid['type'] == 6) & (xid['mode'] == 1) &
                                (xid['Dist(arcsec)'] <= max_radius.to(u.arcsec).value))[0]

            if silent == False:
                print '## Finding nearby stars for '+ID[ii]+'. '+\
                    str(len(good))+' found.'
                if len(good) == 0: print '## Skipping ahead.'

            # Mod on 02/01/2017 to handle python crash when [good] is empty
            if len(good) > 0:
                xid = xid[good]

                # Sort by distance and then brightness
                #xid.sort(['Dist(arcsec)',mag_filt])

                # Sort by brightness and then distance | Mod on 08/01/2017
                xid.sort([mag_filt,'Dist(arcsec)'])

                if silent == False:
                    print '### Writing: ', out_table_file
                    asc.write(xid, out_table_file, overwrite=True,
                              format='fixed_width_two_line')
            else: skip0 = 1
        #endelse

        if len(xid) > 0 and skip0 == 0:
            # Fix bug so c1 is sorted consistently with [xid] | + on 03/01/2017
            # Mod on 24/01/2017 for fk5
            # Moved lower on 26/01/2017
            c1 = coords.SkyCoord(xid['ra'], xid['dec'], 'fk5', unit=(u.deg))

            # Moved up on 28/01/2017
            c1_new, c1_2000 = None, None # Mod on 30/01/017 to handle bug with pmfix=False
            if pmfix == True:
                if ii in m_idx: c1_new, c1_2000 = m_c0[ii], m_c0_2000[ii]
                if ii in u_idx: c1_new, c1_2000 = u_c0[ii], u_c0_2000[ii]
                if ii in s2_idx: c1_new, c1_2000 = s2_c0[ii], s2_c0_2000[ii] # + on 29/01/2017

            if catalog == 'SDSS':
                # + on 9/01/2017
                mag_str, mag_table = sdss_mag_str(xid, TWOMASS=True,
                                                  runall=runall)

            # + on 10/10/2017
            if catalog == '2MASS':
                mag_str = ['JHK=None'] * len(xid)
                for mm in range(len(xid)):
                    tab0 = xid[mm]
                    mag_str[mm] = 'J=%5.2f  H=%5.2f  K=%5.2f' % \
                                  (tab0['j_m'],tab0['h_m'],tab0['k_m'])

                mag_table = xid

            # Always initialize with the proper motion from multi-SDSS and
            # 2MASS | + on 29/01/2017
            # Move out of SDSS if statement, so that proper motion is included
            # for both SDSS and 2MASS | Mod on 10/10/2017
            if pmfix == True:
                mag_table['pRA'][0]    = t_s2['pmRA'][ii]
                mag_table['pDec'][0]   = t_s2['pmDE'][ii]
                mag_table['e_pRA'][0]  = t_s2['e_pmRA'][ii]
                mag_table['e_pDec'][0] = t_s2['e_pmDE'][ii]

            if c1_new != None:
                if ii in m_idx:
                    mag_table['pRA'][0]      = t_movers['pmRA'][ii]
                    mag_table['pDec'][0]     = t_movers['pmDE'][ii]
                    mag_table['e_pRA'][0]    = t_movers['e_pmRA'][ii]
                    mag_table['e_pDec'][0]   = t_movers['e_pmDE'][ii]
                    mag_table['p_source'][0] = 'MoVeRS'
                if ii in u_idx:
                    mag_table['pRA'][0]      = t_ucac['pmRA'][ii]
                    mag_table['pDec'][0]     = t_ucac['pmDE'][ii]
                    mag_table['e_pRA'][0]    = t_ucac['e_pmRA'][ii]
                    mag_table['e_pDec'][0]   = t_ucac['e_pmDE'][ii]
                    # Mod on 09/10/2017
                    if UCAC5 == False:
                        mag_table['p_source'][0] = 'UCAC4'
                    else:
                        mag_table['p_source'][0] = 'UCAC5'

            # + on 10/01/2017
            # Mod on 10/10/2017 to work for both SDSS and 2MASS
            name_Col = Column(np.repeat(ID0[ii]+'_off',len(mag_str)),
                              name='ID')
            mag_table.add_column(name_Col, 0)
            print mag_table0
            if len(mag_table0) == 0: #== None:
                mag_table0 = Table(mag_table[0])
            else:
                mag_table0.add_row(mag_table[0])

            # Moved up on 20/01/2017
            if 'SDSS' in image:
                out_fits = finding_chart_fits_path + ID0[ii]+'.SDSS.fits.gz'
                print out_fits
                if not exists(out_fits):
                    t_hdu = get_sdss_images(c0, out_fits,
                                            band=image.replace('SDSS-',''))

                # + on 06/01/2017
                # Mod on 07/01/2017 to check if FITS file exists
                out_image = finding_chart_fits_path + ID0[ii]+'.crop.SDSS.fits'
                if not exists(out_image):
                    montage_reproj.main(c0, fitsfile=out_fits, catalog=catalog,
                                        out_image=out_image)

            # Moved up on 20/01/2017 | + on 19/01/2017
            if '2MASS' in image:
                out_fits = finding_chart_fits_path + ID0[ii]+'.2MASS.fits.gz'
                if not exists(out_fits): # Mod on 24/01/2017 to over override
                    t_hdu = get_2mass_images(c0, out_fits,
                                             band=image.replace('2MASS-',''))
                out_image = out_fits

            print '### out_fits : ', out_fits

            out_pdf = finding_chart_path + ID0[ii]+'_'+catalog+'_'+image+'.pdf'
            if catalog == 'SDSS' and ('SDSS' in image):
                out_pdf = finding_chart_path + ID0[ii]+'.SDSS.pdf'

            if catalog == '2MASS' and ('2MASS' in image):
                out_pdf  = finding_chart_path + ID0[ii]+'.2MASS.pdf'

            # + on 26/01/2017
            if pmfix == True: out_pdf = out_pdf.replace('.pdf', '.PMfix.pdf')

            print '### out_pdf : ', out_pdf

            # Mod on 02/01/2017 and 22/01/2017 for inputs
            #if catalog == 'SDSS':
            plot_finding_chart(out_image, ID0[ii], band0, c0, c1, mag_str,
                               mag_table, slitlength=slitlength,
                               catalog=catalog, image=image, out_pdf=out_pdf,
                               pmfix=pmfix, c1_new=c1_new, c1_2000=c1_2000,
                               epoch=epoch, MMT=MMT)
        #endif
    #endfor

    # Output table for pmfix case | Mod on 29/01/2017
    if pmfix == False:
        out_mag_table = out_path + 'Alignment_Stars.txt'
    else:
        out_mag_table = alignment_file.replace('Stars.','Stars.PMfix.')

    # + on 10/01/2017. Mod on 29/01/2017
    if catalog == 'SDSS':
        if silent == False: print '### Writing : ', out_mag_table
        mag_table0.write(out_mag_table, format='ascii.fixed_width_two_line',
                         overwrite=True)

    if silent == False:
        print '### End find_nearby_bright_star.main | '+systime()
#enddef

def zcalbase_gal_gemini():
    '''
    Function to run find_nearby_bright_star.main() but for Gemini-N/GNIRS
    targets

    Parameters
    ----------
    None
          
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 23 December 2016
    Modified by Chun Ly, 24 December 2016
     - Added call to main() using 2MASS selection
     - Change max_radius from 5 arcmin to 120 arcsec.
       Slit length of GNIRS is 99 arcsec. Will do offset star if too far away
    Modified by Chun Ly, 2 January 2017
     - Two separate paths for finding chart FITS and PDF
    Modified by Chun Ly, 8 January 2017
     - Use pdfmerge package to merge PDF files
    Modified by Chun Ly, 8 January 2017
     - Run 2MASS before SDSS
    Modified by Chun Ly, 10 January 2017
     - Write photometric ASCII catalog for 2017A GNIRS targets
    Modified by Chun Ly, 20 January 2017
     - Handle making 2MASS finding chart with SDSS coordinates
    Modified by Chun Ly, 25 January 2017
     - Run main() with proper motion information incorporated for finding chart
    '''

    import pdfmerge

    path0                   = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    infile                  = path0 + 'targets.txt'
    out_path                = path0 + 'Alignment_Stars/'
    finding_chart_path      = path0 + 'Finding_Charts/'

    # + on 02/01/2017
    finding_chart_fits_path = '/Users/cly/data/Observing/Gemini/Finding_Charts/'

    max_radius = 95. * u.arcsec #120 * u.arcsec # Mod on 08/01/2017

    slitlength = 99 * u.arcsec 

    # Moved up on 20/01/2017
    infile2 = path0 + 'targets.2017a.txt'
    print '### Reading : ', infile2
    data2   = asc.read(infile2, format='commented_header')

    # Select alignment stars based on 2MASS
    # + on 24/12/2016
    # Moved up on 09/01/2017
    # main(infile, out_path, finding_chart_path, finding_chart_fits_path,
    #     max_radius=max_radius, mag_limit=17.0, catalog='2MASS', mag_filt='j_m')

    do_step1 = 1
    if do_step1:
        # Select alignment stars based on SDSS
        main(infile, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             slitlength=slitlength, runall=True) #runall=False)

        # Merge PDF finding chart files for 2017A targets | + on 08/01/2017
        files = [finding_chart_path+a.replace('*','')+'.SDSS.pdf' for
                 a in data2['ID']]
        #print files
        out_pdf_2017a = finding_chart_path+\
                        'GNIRS_2017A_Targets_SDSS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017a)

    # + on 20/01/2017
    do_step2 = 1
    if do_step2:
        # Generate 2MASS finding chart with SDSS catalog | + on 20/01/2017
        main(infile, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='2MASS-H', slitlength=slitlength, runall=False)

        # Merge PDF finding chart files for 2017A targets
        files = [finding_chart_path+a.replace('*','')+'_SDSS_2MASS-H.pdf' for
                 a in data2['ID']]
        out_pdf_2017a = finding_chart_path+\
                        'GNIRS_2017A_Targets_2MASS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017a)

    # Write Alignment Star photometric catalog
    # + on 10/01/2017
    out_mag_table = out_path + 'Alignment_Stars.txt'
    print '### Reading : ', out_mag_table
    mag_table0 = asc.read(out_mag_table)

    ID1 = [str.replace('_off','') for str in mag_table0['ID']]
    ID2 = [str.replace('*','') for str in data2['ID']]

    idx2, idx1 = match_nosort_str(ID2, ID1)
    out_mag_table2 = out_path + 'Alignment_Stars.2017a.txt'
    print '### Writing : ', out_mag_table2
    asc.write(mag_table0[idx1], out_mag_table2, format='fixed_width_two_line',
              overwrite=True)

    # Run through but use proper motion from UCAC4 or MoVeRS | + on 25/01/2017
    do_step3 = 1
    if do_step3:
        # Generate 2MASS finding chart with SDSS catalog | + on 20/01/2017
        pm_out_pdf = path0 + 'sdss_2mass_proper_motion.pdf'
        main(infile, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='2MASS-H', slitlength=slitlength, runall=False,
             alignment_file=out_mag_table2, pmfix=True, epoch=2017.25,
             pm_out_pdf=pm_out_pdf)

        # Merge PDF finding chart files for 2017A targets
        files = [finding_chart_path+a.replace('*','')+'_SDSS_2MASS-H.PMfix.pdf' for
                 a in data2['ID']]
        out_pdf_2017a = finding_chart_path+\
                        'GNIRS_2017A_Targets_2MASS_FindingCharts.PMfix.pdf'
        pdfmerge.merge(files, out_pdf_2017a)

        ## Get updated table with proper motion | + on 29/01/2017
        #out_mag_table = out_path + 'Alignment_Stars.PMfix.txt'
        #print '### Reading : ', out_mag_table
        #mag_table0 = asc.read(out_mag_table)
        #
        #out_mag_table2 = out_path + 'Alignment_Stars.PMfix.2017a.txt'
        #print '### Writing : ', out_mag_table2
        #asc.write(mag_table0[idx1], out_mag_table2, overwrite=True,
        #          format='fixed_width_two_line')
#enddef

def zcalbase_gal_gemini_2017b():
    '''
    Function to run find_nearby_bright_star.main() but for 2017B Gemini-N/GNIRS
    targets

    Parameters
    ----------
    None

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 7 July 2017
     - Started as a copy of zcalbase_gal_gemini()
     - Minor fixes. Read in alignment table
    Modified by Chun Ly, 8 July 2017
     - Adopt 2.5 sigma constraint on PM reliability
     - Change epoch for 2017B observation planning: Using Nov 3 2017 as Field4
       transits at midnight
    '''

    import pdfmerge

    path0                   = '/Users/cly/Dropbox/Observing/2017B/Gemini/'
    out_path                = path0 + 'Alignment_Stars/'
    finding_chart_path      = path0 + 'Finding_Charts/'

    finding_chart_fits_path = '/Users/cly/data/Observing/Gemini/Finding_Charts/'

    max_radius = 95. * u.arcsec
    slitlength = 99 * u.arcsec

    infile2 = path0 + 'targets.2017b.txt'
    print '### Reading : ', infile2
    data2   = asc.read(infile2, format='commented_header')

    do_step1 = 1
    if do_step1:
        # Select alignment stars based on SDSS
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             slitlength=slitlength, runall=True) #runall=False)

        # Merge PDF finding chart files for 2017A targets
        files = [finding_chart_path+a.replace('*','')+'.SDSS.pdf' for
                 a in data2['ID']]

        out_pdf_2017b = finding_chart_path+\
                        'GNIRS_2017B_Targets_SDSS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

    do_step2 = 1
    if do_step2:
        # Generate 2MASS finding chart with SDSS catalog
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='2MASS-H', slitlength=slitlength, runall=False)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'_SDSS_2MASS-H.pdf' for
                 a in data2['ID']]
        out_pdf_2017b = finding_chart_path+\
                        'GNIRS_2017B_Targets_2MASS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

    out_mag_table = out_path + 'Alignment_Stars.txt'
    print '### Reading : ', out_mag_table
    mag_table0 = asc.read(out_mag_table)

    # Run through but use proper motion from UCAC4 or MoVeRS
    do_step3 = 1
    if do_step3:
        # Generate 2MASS finding chart with SDSS catalog
        pm_out_pdf = path0 + 'sdss_2mass_proper_motion.pdf'
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='2MASS-H', slitlength=slitlength, runall=False,
             alignment_file=out_mag_table, pmfix=True, sig_min=2.5,
             epoch=2017.84, pm_out_pdf=pm_out_pdf)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'_SDSS_2MASS-H.PMfix.pdf' for
                 a in data2['ID']]
        out_pdf_2017b = finding_chart_path+\
                        'GNIRS_2017B_Targets_2MASS_FindingCharts.PMfix.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

#enddef


def zcalbase_gal_mmt_2017b():
    '''
    Function to run find_nearby_bright_star.main() but for 2017B MMT/MMIRS
    targets

    Parameters
    ----------
    None

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 30 July 2017
     - Started as a copy of zcalbase_gal_gemini_2017b()
    Modified by Chun Ly, 8 October 2017
     - No longer generate 2MASS finding charts with proper motion fix
     - Fix typo bug
     - Handle MMT case in main()
    Modified by Chun Ly, 9 October 2017
     - Use UCAC5 proper motion
     - Rename output sdss_2mass proper motion file
     - Use proper length for longslit (~7 arcmin)
    Modified by Chun Ly, 10 October 2017
     - do_step4: Run with 2MASS catalog input without proper motion
    '''

    import pdfmerge

    path0                   = '/Users/cly/Dropbox/Observing/2017B/MMT/'
    out_path                = path0 + 'Alignment_Stars/'
    finding_chart_path      = path0 + 'Finding_Charts/'

    finding_chart_fits_path = '/Users/cly/data/Observing/MMIRS/Finding_Charts/'

    # Mod on 09/10/2017
    max_radius = 120. * u.arcsec
    slitlength = 2048 * 0.201 * u.arcsec #4 * 60 * u.arcsec

    infile2 = path0 + 'targets.2017b.txt'
    print '### Reading : ', infile2
    data2   = asc.read(infile2, format='commented_header')

    do_step1 = 1
    if do_step1:
        # Select alignment stars based on SDSS
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             slitlength=slitlength, runall=True, MMT=True) #runall=False)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'.SDSS.pdf' for
                 a in data2['ID']]

        out_pdf_2017b = finding_chart_path+\
                        'MMIRS_2017B_Targets_SDSS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

    do_step2 = 1
    if do_step2:
        # Generate 2MASS finding chart with SDSS catalog
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='2MASS-H', slitlength=slitlength, runall=False, MMT=True)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'_SDSS_2MASS-H.pdf' for
                 a in data2['ID']]
        out_pdf_2017b = finding_chart_path+\
                        'MMIRS_2017B_Targets_2MASS_FindingCharts.bright.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

    out_mag_table = out_path + 'Alignment_Stars.txt'
    print '### Reading : ', out_mag_table
    mag_table0 = asc.read(out_mag_table)

    # Run through but use proper motion from UCAC4 or MoVeRS
    do_step3 = 1
    if do_step3:
        # Generate SDSS finding chart with SDSS catalog
        pm_out_pdf = path0 + 'sdss_2mass_proper_motion.UCAC5.pdf'
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
             image='SDSS', slitlength=slitlength, runall=False,
             alignment_file=out_mag_table, pmfix=True, sig_min=2.5,
             epoch=2017.84, pm_out_pdf=pm_out_pdf, MMT=True, UCAC5=True)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'.SDSS.PMfix.pdf' for
                 a in data2['ID']]
        out_pdf_2017b = finding_chart_path+\
                        'MMIRS_2017B_Targets_SDSS_FindingCharts.PMfix.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

    do_step4 = 1
    if do_step4:
        # Generate SDSS finding chart with 2MASS catalog with PM
        pm_out_pdf = path0 + 'sdss_2mass_proper_motion_2MASS.UCAC5.pdf'
        out_mag_table = out_path + 'Alignment_Stars_2MASS.txt'
        main(infile2, out_path, finding_chart_path, finding_chart_fits_path,
             max_radius=max_radius, mag_limit=19.0, mag_filt='j_m',
             catalog='2MASS', image='SDSS', slitlength=slitlength, runall=False,
             alignment_file=out_mag_table, pmfix=False, sig_min=2.5,
             epoch=2017.84, pm_out_pdf=pm_out_pdf, MMT=True, UCAC5=True)

        # Merge PDF finding chart files for 2017B targets
        files = [finding_chart_path+a.replace('*','')+'_2MASS_SDSS.pdf' for
                 a in data2['ID']]
        out_pdf_2017b = finding_chart_path+\
                        'MMIRS_2017B_Targets_2MASS-SDSS_FindingCharts.pdf'
        pdfmerge.merge(files, out_pdf_2017b)

#enddef
