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
SDSS_fld      = ['ra','dec','objid','run','rerun','camcol','field','obj',
                 'type','mode']
                 
SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

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

def main(infile, out_path, finding_chart_path, max_radius=60*u.arcsec,
         mag_limit=20.0, mag_filt='modelMag_i', catalog='SDSS',
         format='commented_header', silent=False, verbose=True):

    '''
    Main function to find nearby star

    Parameters
    ----------
    infile : string
      Filename for ASCII file that contains ID, RA, and DEC
      RA and DEC must be provided in sexigesimal format.
      The columns should be referred to as 'ID', 'RA', and 'DEC'

    out_path : string
      Full path to output ASCII tables.  Must end with a '/'

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
    
    for ii in range(3): #n_sources):
        c0 = coords.SkyCoord(ra=RA[ii], dec=DEC[ii], unit=(u.hour, u.degree))
        if catalog == 'SDSS':
            xid = SDSS.query_region(c0, max_radius, data_release=12,
                                    photoobj_fields=SDSS_phot_fld)

            # Keep stars only (type == 6)
            # http://www.sdss.org/dr12/algorithms/classify/#photo_class
            # Keep primary target to deal with duplicate entries | 24/12/2016
            good = np.where((xid[mag_filt] <= mag_limit) &
                            (xid['type'] == 6) & (xid['mode'] == 1))[0]
            # Moved up on 24/12/2016
            out_table_file = out_path + ID0[ii] + '.SDSS.nearby.txt'
            
        # + on 24/12/2016
        if catalog == '2MASS':
            xid = IRSA.query_region(c0, catalog='fp_psc', radius=max_radius)
            good = np.where(xid[mag_filt] <= mag_limit)[0]
            out_table_file = out_path + ID0[ii] + '.2MASS.nearby.txt'

        xid = xid[good]

        if silent == False:
            print '## Finding nearby stars for '+ID[ii]+'. '+\
                str(len(good))+' found.'

        # Get distance from target
        c1 = coords.SkyCoord(xid['ra'], xid['dec'], unit=(u.degree, u.degree))
        sep = c0.separation(c1).to(u.arcsec).value
        col1 = Column(sep, name='Dist(arcsec)')
        xid.add_column(col1)

        # Sort by distance and then brightness
        xid.sort(['Dist(arcsec)',mag_filt])
        
        if silent == False:
            print '### Writing: ', out_table_file
        asc.write(xid, out_table_file, format='fixed_width_two_line')

        if catalog == 'SDSS':
            out_fits = finding_chart_path + ID0[ii]+'.SDSS.fits.gz'
            print out_fits
            if not exists(out_fits):
                t_hdu = get_sdss_images(c0, out_fits,
                                        band=mag_filt.replace('modelMag_',''))
            else:
                t_hdu = fits.open(out_fits)
    #endfor

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
    '''

    path0              = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    infile             = path0 + 'targets.txt'
    out_path           = path0 + 'Alignment_Stars/'
    finding_chart_path = path0 + 'Finding_Charts/'

    # Select alignment stars based on SDSS
    main(infile, out_path, finding_chart_path, max_radius=5*u.arcmin,
         mag_limit=19.0, catalog='SDSS')

    # Select alignment stars based on 2MASS
    # + on 24/12/2016
    main(infile, out_path, finding_chart_path, max_radius=5*u.arcmin,
         mag_limit=17.0, catalog='2MASS', mag_filt='j_m')

