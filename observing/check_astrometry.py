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

from astroquery.sdss import SDSS # + on 12/01/2017

# Box region around target's coordinates to perform astrometry check
# + on 12/01/2017
box_size = (5 * u.arcmin).to(u.arcsec).value

SDSS_fld      = ['ra','dec','objid','run','rerun','camcol','field','obj',
                 'type','mode']

SDSS_phot_fld = SDSS_fld + ['modelMag_u', 'modelMagErr_u', 'modelMag_g',
                            'modelMagErr_g', 'modelMag_r', 'modelMagErr_r',
                            'modelMag_i', 'modelMagErr_i', 'modelMag_z',
                            'modelMagErr_z']

def main(c0, label, infile1=None, data1=None, infile2=None, cat2_survey='SDSS',
         columns=['ALPHA_J2000', 'DELTA_J2000'], out_pdf=None, silent=False,
         verbose=True):

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
    if infile1 == None and data1 == None:
        print '''
        '## Error. Must provide infile1 or data1'
        '## Exiting'
        '''
        return
    #endif

    rad0 = 1; src_rad0 = rad0 * u.arcsec

    # + on 12/01/2017
    if infile1 != None:
        if silent == False: print '### Reading : ', infile1
        if '.fits' in infile1:
            hdu1  = fits.open(infile1)
            data1 = hdu1[1].data
            cmd1a = 'RA  = data1.'+columns[0]
            cmd1b = 'DEC = data1.'+columns[1]
            exec(cmd1a); exec(cmd1b)
        else:
            data1 = asc.read(infile1)
            RA  = data1[columns[0]]
            DEC = data1[columns[1]]

        RA_off  = (RA1 - c0.ra.deg) * 3600.0 * np.cos(np.radians(c0.dec.deg))
        DEC_off = (DEC1 - c0.dec.deg) * 3600.0
        in_box  = np.where((np.absolute(RA_off) <= box_size) &
                           (np.absolute(DEC_off) <= box_size))[0]
        data1   = data1[in_box]
    else:
        if silent == False: print '## [data1] defined'

    # Get coordinates again | + on 12/01/2017
    cmd1a = 'RA  = data1.'+columns[0]
    cmd1b = 'DEC = data1.'+columns[1]
    exec(cmd1a); exec(cmd1b)

    c1_arr = coords.SkyCoord(ra=RA, dec=DEC, unit=(u.deg,u.deg)) # + on 12/01/2017

    # Get reference catalog | + on 12/01/2017
    if infile2 == None:
        if cat2_survey == 'SDSS':
            s_rad0 = box_size*u.arcsec*np.sqrt(2) # radius in arcsec
            ref_cat0 = SDSS.query_region(c0, radius=s_rad0, data_release=12,
                                    photoobj_fields=SDSS_phot_fld)
    else:
        if silent == False: print '### Reading : ', infile2
        ref_cat0 = fits.getdata(infile2)

    if verbose == True:
        print '## '+cat2_survey+' size: [ref_cat0] ', len(ref_cat0)

    # + on 12/01/2017
    if cat2_survey == 'SDSS':
        RA_ref  = ref_cat0['ra']
        DEC_ref = ref_cat0['dec']
        c1_ref  = coords.SkyCoord(ra=RA_ref, dec=DEC_ref, unit=(u.deg,u.deg))

        #idx, d2d, d3d = c1_arr.match_to_catalog_sky(c1_ref)
        #print idx; print d2d; print d3d
        idx_arr, idx_ref, d2d, d3d = c1_ref.search_around_sky(c1_arr, src_rad0)
        print len(idx_ref), len(idx_arr), len(d2d), len(d3d)

        A_RA, A_DEC = RA[idx_arr], DEC[idx_arr]
        R_RA, R_DEC = RA_ref[idx_ref], DEC_ref[idx_ref]
        diff_RA  = (A_RA - R_RA)   * 3600.0 * np.cos(np.radians(A_DEC))
        diff_DEC = (A_DEC - R_DEC) * 3600.0

    # Generate plot from cross-match results | + on 12/01/2017
    fig, ax = plt.subplots()
    ax.scatter(diff_RA, diff_DEC, marker='o', c='none', s=50,
               edgecolors='black')
    ax.set_xlim(rad0*np.array([-1,1]))
    ax.set_ylim(rad0*np.array([-1,1]))
    ax.minorticks_on()
    ax.set_xlabel('RA: Primary - '+cat2_survey+' [arcsec]', fontsize='16')
    ax.set_ylabel('DEC: Primary - '+cat2_survey+' [arcsec]', fontsize='16')
    ax.annotate(label, [0.025,0.975], xycoords='axes fraction', ha='left',
                va='top', fontsize='large', weight='bold')
    fig.set_size_inches(8,8)
    fig.tight_layout()

    if infile1 != None and out_pdf == None:
        out_pdf = infile1.replace('.fits.gz','.pdf').replace('.dat','.pdf')

    if silent == False: print '### Writing :', out_pdf
    fig.savefig(out_pdf)

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

    path0      = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    cat_path0  = '/Users/cly/data/SDF/NBcat/SEx/'
    astro_path = path0 + 'Astrometry/'

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

        out_pdf = astro_path + ID[ii]+'.astrometry.pdf'
        main(t_c, ID[ii], data1=box_data1, columns=['ALPHA_J2000','DELTA_J2000'],
             out_pdf=out_pdf, verbose=True)

    #endfor
    if silent == False: print '### End check_astrometry.SDF | '+systime()
#enddef
