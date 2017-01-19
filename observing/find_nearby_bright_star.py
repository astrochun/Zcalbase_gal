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

import aplpy # + on 24/12/2016

from PyMontage.scripts import montage_reproj # + on 02/01/2017

# For SDSS only
SDSS_fld      = ['ra','dec','objid','run','rerun','camcol','field','obj',
                 'type','mode']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def get_PA(c0, c1, slitlength=99*u.arcsec, silent=True, verbose=False):
    '''
    Function to determine the PA on the sky and central coordinates given
    two WCS coordinates

    Parameters
    ----------
    c0 : `astropy.coordinates` object
      Central coordinate of target

    c1 : `astropy.coordinates` object
      Central coordinate of nearby bright star

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
    '''
    
    if silent == False:
        print '### Begin find_nearby_bright_star.get_PA | '+systime()

    PA  = c0.position_angle(c1).degree # + => East of North
    PAr = np.radians(PA)

    # Get central coordinate
    ra_avg  = np.average([c0.ra.value, c1.ra.value])
    dec_avg = np.average([c0.dec.value, c1.dec.value])

    c_ctr = coords.SkyCoord(ra=ra_avg, dec=dec_avg, unit=(u.degree,u.degree))

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
    '''

    # Give coordinate offsets from reference [c_ref] to target [c0]
    dra0, ddec0 = c_ref.spherical_offsets_to(c0)
    dra0  = dra0.to(u.arcsec).value
    ddec0 = ddec0.to(u.arcsec).value
    dist0 = np.sqrt(dra0**2+ddec0**2)
    return dra0, ddec0, dist0
#enddef

def plot_finding_chart(fitsfile, t_ID, band0, c0, c1, mag_str, out_pdf=None,
                       slitlength=99*u.arcsec, catalog='SDSS', silent=True,
                       verbose=False):
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

    out_pdf : string
      Output PDF filename.
      Default: based on [fitsfile], replacing '.fits.gz' or '.fits' with '.pdf'

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
    '''

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
    PA, c_ctr, longslit_list = get_PA(c0, c1[0], slitlength=slitlength)

    gc.show_lines(longslit_list, layer='slit', color='black', linewidth=0.5)
    #gc.show_rectangles([c_ctr.ra.value], [c_ctr.dec.value], 1/3600.0, 99/3600.0)

    # Fix bug. marker='+' won't work with facecolor='none'
    gc.show_markers([c0.ra.value], [c0.dec.value], layer='primary',
                    edgecolor='red', facecolor='red', marker='+', s=25)

    # + on 02/01/2017
    gc.show_markers([c1.ra.value], [c1.dec.value], layer='secondary',
                    edgecolor='blue', facecolor='none', marker='o', s=25,
                    linewidth=0.5)
    
    # Label things in lower left text | + on 03/01/2017
    str_c_t  = c0.to_string('hmsdms').split(' ')
    str_c    = c_ctr.to_string('hmsdms').split(' ')
    str_c_bt = c1[0].to_string('hmsdms').split(' ')
    dra, ddec, dist    = get_offsets(c1[0], c0)    # Mod on 10/01/2017
    dra2, ddec2, dist2 = get_offsets(c1[0], c_ctr) # Mod on 10/01/2017

    # Mod on 09/01/2017, 10/01/2017
    bt_txt  = 'Target: RA='+str_c_t[0]+', Dec='+str_c_t[1]+'\n\n'
    bt_txt += 'Slit Center: RA='+str_c[0]+', Dec='+str_c[1]+'\n'
    bt_txt += ('Slit PA = %7.3f' % PA) + ' deg\n'
    bt_txt += 'Offsets : (%+.3f", %+.3f");  %.2f"\n\n' % (dra2, ddec2, dist2)
    bt_txt += 'Offset Star: RA='+str_c_bt[0]+', Dec='+str_c_bt[1]+'\n'
    bt_txt += 'Offsets : (%+.3f", %+.3f");  %.2f"\n' % (dra, ddec, dist)
    bt_txt += mag_str[0]

    gc.add_label(0.03, 0.125, bt_txt, color='magenta', relative=True,
                 ha='left', va='bottom', weight='medium', size='small')

    # Label upper left source, catalog, and band | + on 02/01/2017
    lab0 = t_ID+'\n'+catalog+' '+band0
    gc.add_label(0.03, 0.95, lab0, relative=True, ha='left', va='top',
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
    '''

    if silent == False:
        print '### Begin find_nearby_bright_star.get_2mass_images '+systime()

    imgs = SkyView.get_images(coordinates=c0, band='2MASS-'+band,
                              width=5*u.arcmin, height=5*u.arcmin,
                              timeout=180)

    if silent == False: print '### Writing : ', out_fits
    for ff in range(n_frames):
        t_hdu = imgs[0][0]
        if ff == 0:
            t_hdu.writeto(out_fits, clobber=True)
        else:
            fits.append(out_fits, t_hdu.data, t_hdu.header)

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

    if TWOMASS == True:
        t_c0 = coords.SkyCoord(ra=table['ra'], dec=table['dec'],
                               unit=(u.deg, u.deg))

        n_val = n_sources if runall == True else 1
        for ii in range(n_val):
            table_2mass = IRSA.query_region(t_c0[ii], catalog='fp_psc',
                                            radius=1*u.arcsec)
            if len(table_2mass) == 0:
                mag_str[ii] += 'JHK=None'
            else:
                tab0 = table_2mass[0]
                mag_str[ii] += 'J=%5.2f  H=%5.2f  K=%5.2f' % \
                               (tab0['j_m'],tab0['h_m'],tab0['k_m'])

                # + on 10/01/2017
                col_Jmag[ii]  = tab0['j_m']
                col_Hmag[ii]  = tab0['h_m']
                col_Kmag[ii]  = tab0['k_m']
                col_eJmag[ii] = tab0['j_cmsig']
                col_eHmag[ii] = tab0['h_cmsig']
                col_eKmag[ii] = tab0['k_cmsig']
        #endfor
    #endif

    # + on 10/01/2017
    mag_table = Table(mag_table)
    n_cols = len(mag_table.colnames)
    mag_table.add_columns([col_Jmag, col_eJmag, col_Hmag, col_eHmag,
                           col_Kmag, col_eKmag], np.repeat(n_cols-1, 6))
    return mag_str, mag_table
#enddef

def main(infile, out_path, finding_chart_path, finding_chart_fits_path,
         max_radius=60*u.arcsec, mag_limit=20.0, mag_filt='modelMag_i',
         catalog='SDSS', format='commented_header', slitlength=99*u.arcsec,
         runall=True, silent=False, verbose=True):

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
      Either SDSS or '2MASS'. Default: 'SDSS'
    
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
    '''

    if silent == False:
        print '### Begin find_nearby_bright_star.main | '+systime()

    if silent == False: print '### Reading : ', infile
    data0 = asc.read(infile)

    ID  = data0['ID'].data
    RA  = data0['RA'].data
    DEC = data0['DEC'].data
    n_sources = len(ID)

    ID0 = [ss.replace('*','') for ss in ID] # Clean things up 24/12/2016

    if silent == False:
        print '## Search criteria : '
        print '## max_radius(arcsec) : ', max_radius.to(u.arcsec).value
        print '## mag_limit : ', mag_limit
        print '## filter selection : ', mag_filt

    mag_table0 = None # Initialize | + on 10/01/2017
    for ii in range(n_sources):
        c0 = coords.SkyCoord(ra=RA[ii], dec=DEC[ii], unit=(u.hour, u.degree))
        # Moved up on 09/01/2017
        out_table_file = out_path + ID0[ii] + '.'+catalog+'.nearby.txt' 

        if catalog == 'SDSS':
            xid = SDSS.query_region(c0, radius=max_radius, data_release=12,
                                    photoobj_fields=SDSS_phot_fld)

        # + on 24/12/2016 | Moved up on 09/01/2017
        if catalog == '2MASS':
            xid = IRSA.query_region(c0, catalog='fp_psc', radius=max_radius)
            good = np.where(xid[mag_filt] <= mag_limit)[0]

        # Get distance from target
        if len(xid) > 0:
            c1 = coords.SkyCoord(xid['ra'], xid['dec'], unit=(u.deg, u.deg))
            sep = c0.separation(c1).to(u.arcsec).value
            col1 = Column(sep, name='Dist(arcsec)')
            xid.add_column(col1)

        if catalog == 'SDSS':
            # Keep stars only (type == 6)
            # http://www.sdss.org/dr12/algorithms/classify/#photo_class
            # Keep primary target to deal with duplicate entries | 24/12/2016
            # Avoid unusual mag values | Mod on 24/12/2016
            # Fix SDSS.query_region problem with max_radius | Mod on 08/01/2017
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

            # Fix bug so c1 is sorted consistently with [xid] | + on 03/01/2017
            c1 = coords.SkyCoord(xid['ra'], xid['dec'], unit=(u.deg, u.deg))

            if silent == False:
                print '### Writing: ', out_table_file
                asc.write(xid, out_table_file, format='fixed_width_two_line',
                          overwrite=True)
        
            if catalog == 'SDSS':
                # + on 9/01/2017
                mag_str, mag_table = sdss_mag_str(xid, TWOMASS=True, runall=runall)

                # + on 10/01/2017
                name_Col = Column(np.repeat(ID0[ii]+'_off',len(mag_str)), name='ID')
                mag_table.add_column(name_Col, 0)
                if mag_table0 == None:
                    mag_table0 = Table(mag_table[0])
                else:
                    mag_table0.add_row(mag_table[0])

                out_fits = finding_chart_fits_path + ID0[ii]+'.SDSS.fits.gz'
                out_pdf  = finding_chart_path + ID0[ii]+'.SDSS.pdf'
                print out_fits
                band0 = mag_filt.replace('modelMag_','') # + on 02/01/2017
                if not exists(out_fits):
                    t_hdu = get_sdss_images(c0, out_fits, band=band0)
                else:
                    t_hdu = fits.open(out_fits)

                # + on 06/01/2017
                # Mod on 07/01/2017 to check if FITS file exists
                out_image = finding_chart_fits_path + ID0[ii]+'.crop.SDSS.fits'
                if not exists(out_image):
                    montage_reproj.main(c0, fitsfile=out_fits, catalog=catalog,
                                        out_image=out_image)

            # Mod on 02/01/2017 for inputs
            if catalog == 'SDSS':
                plot_finding_chart(out_image, ID0[ii], band0, c0, c1,
                                   mag_str, slitlength=slitlength,
                                   catalog=catalog, out_pdf=out_pdf)
        #endif
    #endfor

    # + on 10/01/2017
    out_mag_table = out_path + 'Alignment_Stars.txt'
    if silent == False:
        print '### Writing : ', out_mag_table
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

    # Select alignment stars based on 2MASS
    # + on 24/12/2016
    # Moved up on 09/01/2017
    #main(infile, out_path, finding_chart_path, finding_chart_fits_path,
    #     max_radius=max_radius, mag_limit=17.0, catalog='2MASS', mag_filt='j_m')

    # Select alignment stars based on SDSS
    main(infile, out_path, finding_chart_path, finding_chart_fits_path,
         max_radius=max_radius, mag_limit=19.0, catalog='SDSS',
         slitlength=slitlength, runall=False)

    # Merge PDF finding chart files for 2017A targets | + on 08/01/2017
    infile2 = path0 + 'targets.2017a.txt'
    print '### Reading : ', infile2
    data2   = asc.read(infile2, format='commented_header')
    files = [finding_chart_path+a.replace('*','')+'.SDSS.pdf' for
             a in data2['ID']]
    #print files
    out_pdf_2017a = finding_chart_path+'GNIRS_2017A_Targets_FindingCharts.bright.pdf'
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
