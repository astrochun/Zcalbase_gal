"""
check_astrometry
================

Compares the WCS for two different catalogs, such as a target sample and a
reference survey (e.g., 2MASS, SDSS) to check relative astrometry

Requirements:
 astroquery (https://astroquery.readthedocs.io/en/latest/)
 pdfmerge ('pip install pdfmerge' will install it)

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
from astroquery.irsa import Irsa as IRSA # + on 14/01/2017

# Moved up on 12/01/2017
path0      = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
astro_path = path0 + 'Astrometry/'

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
    Modified by Chun Ly, 13 January 2017
     - Output average, median, and sigma of offsets
    Modified by Chun Ly, 14 January 2017
     - Add option for cat2_survey == '2MASS'
     - Add option for cat2_survey == 'WISE'
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

    limits0 = rad0 * np.array([-1,1]) # + on 12/01/2017

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
        # + on 14/01/2017
        if cat2_survey == '2MASS':
            size0 = box_size*u.arcsec # size of edge in arcsec
            ref_cat0 = IRSA.query_region(c0, catalog='fp_psc', spatial='Box',
                                         width=size0)

        # + on 14/01/2017
        if cat2_survey == 'WISE':
            size0 = box_size*u.arcsec # size of edge in arcsec
            ref_cat0 = IRSA.query_region(c0, catalog='allwise_p3as_psd',
                                         spatial='Box', width=size0)

    else:
        if silent == False: print '### Reading : ', infile2
        ref_cat0 = fits.getdata(infile2)

    if verbose == True:
        print '## '+cat2_survey+' size : [ref_cat0] ', len(ref_cat0)

    # + on 12/01/2017
    if cat2_survey == 'SDSS':
        RA_ref  = ref_cat0['ra']
        DEC_ref = ref_cat0['dec']
        c1_ref  = coords.SkyCoord(ra=RA_ref, dec=DEC_ref, unit=(u.deg,u.deg))

    # + on 14/01/2017
    if cat2_survey == '2MASS' or cat2_survey == 'WISE':
        c1_ref  = coords.SkyCoord(ra=ref_cat0['clon'], dec=ref_cat0['clat'],
                                  unit=(u.hour,u.deg))
        RA_ref  = c1_ref.ra.degree
        DEC_ref = c1_ref.dec.degree

    #idx, d2d, d3d = c1_arr.match_to_catalog_sky(c1_ref)
    #print idx; print d2d; print d3d
    idx_arr, idx_ref, d2d, d3d = c1_ref.search_around_sky(c1_arr, src_rad0)
    if silent == False: print '## Crossmatch sample size : ', len(idx_ref)

    A_RA, A_DEC = RA[idx_arr], DEC[idx_arr]
    R_RA, R_DEC = RA_ref[idx_ref], DEC_ref[idx_ref]
    diff_RA  = (A_RA - R_RA)   * 3600.0 * np.cos(np.radians(A_DEC))
    diff_DEC = (A_DEC - R_DEC) * 3600.0

    a_RA, m_RA, s_RA    = np.average(diff_RA), np.median(diff_RA), \
                          np.std(diff_RA)
    a_DEC, m_DEC, s_DEC = np.average(diff_DEC), np.median(diff_DEC), \
                          np.std(diff_DEC)

    # Generate plot from cross-match results | + on 12/01/2017
    fig, ax = plt.subplots()

    str1 = r'RA:  avg=%5.2f med=%5.2f $\sigma$=%5.2f' % (a_RA, m_RA, s_RA)
    str2 = r'DEC: avg=%5.2f med=%5.2f $\sigma$=%5.2f' % (a_DEC, m_DEC, s_DEC)
    str0 = 'N = '+str(len(idx_ref))+'\n'+str1+'\n'+str2

    ax.annotate(str0, [0.975,0.025], xycoords='axes fraction', ha='right',
                va='bottom', fontsize=14)

    ax.plot(limits0, np.repeat(0,2), 'r--')
    ax.plot(np.repeat(0,2), limits0, 'r--')
    ax.scatter(diff_RA, diff_DEC, marker='+', c='k', s=50)
    #ax.scatter(diff_RA, diff_DEC, marker='o', c='none', s=50,
    #edgecolors='black')

    ax.set_xlim(limits0)
    ax.set_ylim(limits0)
    ax.minorticks_on()
    ax.set_xlabel('RA: Primary - '+cat2_survey+' [arcsec]', fontsize='16')
    ax.set_ylabel('DEC: Primary - '+cat2_survey+' [arcsec]', fontsize='16')
    ax.annotate(label, [0.025,0.975], xycoords='axes fraction', ha='left',
                va='top', fontsize='x-large', weight='bold')
    fig.set_size_inches(8,8)
    fig.tight_layout()

    if infile1 != None and out_pdf == None:
        out_pdf = infile1.replace('.fits.gz','.pdf').replace('.dat','.pdf')

    if silent == False: print '### Writing :', out_pdf
    fig.savefig(out_pdf)

    if silent == False: print '### End check_astrometry.main | '+systime()

    # + on 13/01/2017
    return a_RA, m_RA, s_RA, a_DEC, m_DEC, s_DEC
#enddef

def run_main(sample, cat2_survey='SDSS', silent=False, verbose=True):
    '''
    Cross-match either SDF NB catalogs or DEEP2 catalogs against SDSS to get
    astrometry around each [OIII]4363 targets from Ly et al. (2016), ApJS,
    226, 5 and Ly et al. (2015), ApJ, 805, 45

    Parameters
    ----------
    sample : string
      Either 'SDF' or 'DEEP2'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 12 January 2017
    Modified by Chun Ly, 13 January 2017
     - Get average, median, sigma from offsets and create astropy.Table
    Modified by Chun Ly, 14 January 2017
     - Include cat2_survey keyword option.
     - Pass cat2_survey to main()
     - Write ASCII table containing astrometric offsets
    '''

    if silent == False: print '### Begin check_astrometry.run_main | '+systime()

    if sample == 'SDF':
        cat_path0 = '/Users/cly/data/SDF/NBcat/SEx/catalogs/'
        infile0   = path0 + 'targets_sdf_astro.2017a.txt'
    if sample == 'DEEP2':
        cat_path0 = '/Users/cly/data/DEEP2/DR4/phot_cat/'
        infile0   = path0 + 'targets_deep2_astro.2017a.txt'

    if silent == False: print '### Reading : ', infile0
    data0   = asc.read(infile0, format='commented_header')
    print data0

    n_sources = len(data0)

    ID  = data0['ID']
    RA  = data0['RA']
    DEC = data0['DEC']
    c0  = coords.SkyCoord(ra=RA, dec=DEC, unit=(u.hour, u.deg))

    # Mod on 12/01/2017
    if sample == 'SDF':
        cat_files = [cat_path0+'sdf_pub2_'+f0+'.cat.fits.gz' for
                     f0 in data0['filt']]
        #cat_files = [cat_path0+f0+'emitters_allphot.fits' for f0 in data0['filt']]
    if sample == 'DEEP2':
        cat_files = [cat_path0+f0+'_ext.fits.gz' for f0 in data0['field']]
    print cat_files

    # + on 13/01/2017
    a_RA0  = []; m_RA0  = []; s_RA0  = []
    a_DEC0 = []; m_DEC0 = []; s_DEC0 = []
    for ii in range(n_sources):
        if silent == False: print '### Reading : ', cat_files[ii]
        data1 = fits.getdata(cat_files[ii])
        if sample == 'SDF':
            RA1   = data1.ALPHA_J2000
            DEC1  = data1.DELTA_J2000
            columns = ['ALPHA_J2000','DELTA_J2000']
        if sample == 'DEEP2':
            RA1   = data1.RA_SDSS
            DEC1  = data1.DEC_SDSS
            columns = ['RA_SDSS','DEC_SDSS']

        t_c = c0[ii]
        RA_off  = (RA1 - t_c.ra.deg) * 3600.0 * np.cos(np.radians(t_c.dec.deg))
        DEC_off = (DEC1 - t_c.dec.deg) * 3600.0
        in_box = np.where((np.absolute(RA_off) <= box_size) &
                          (np.absolute(DEC_off) <= box_size))[0]
        box_data1 = data1[in_box]

        if silent == False: print '## in_box : ', len(in_box)

        out_pdf = astro_path + ID[ii]+'.'+cat2_survey+'.astrometry.pdf'

        # Mod on 14/01/2017 to include cat2_survey
        val0 = main(t_c, ID[ii], data1=box_data1, cat2_survey=cat2_survey,
                    columns=columns, out_pdf=out_pdf, verbose=True)

        # + on 13/01/2017
        a_RA0.append(val0[0]); m_RA0.append(val0[1]); s_RA0.append(val0[2])
        a_DEC0.append(val0[3]); m_DEC0.append(val0[4]); s_DEC0.append(val0[5])
    #endfor

    # + on 13/01/2017
    names0 = ('ID','avg_RA','med_RA','sig_RA','avg_DEC','med_DEC','sig_DEC')
    tab0 = Table([ID, a_RA0,m_RA0,s_RA0, a_DEC0,m_DEC0,s_DEC0], names=names0)
    print tab0

    # + on 14/01/2017
    out_table_file = astro_path + 'astro_fixes_'+sample+'.'+cat2_survey+'.tbl'
    if silent == False: print '### Writing : ', out_table_file

    tab0.write(out_table_file, format='ascii.fixed_width_two_line')

    if silent == False: print '### End check_astrometry.run_main | '+systime()
#enddef

# Moved up on 14/01/2017
def GNIRS_2017A():
    '''
    Run run_main() for both SDF and DEEP2 to get astrometry plot checks for
    Gemini-N/GNIRS 2017A targets

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
    Modified by Chun Ly, 14 January 2017
     - Specifically set cat2_survey for run_main()
     - Crossmatch against 2MASS option
     - Crossmatch against WISE option
    '''

    import pdfmerge

    ## Do astrometry check against SDSS
    ## Mod on 14/01/2017
    do_SDSS = 1
    if do_SDSS:
        run_main('SDF', cat2_survey='SDSS')
        run_main('DEEP2', cat2_survey='SDSS')

        infile2 = path0 + 'targets.2017a.txt'
        print '### Reading : ', infile2
        data2   = asc.read(infile2)
        files = [astro_path+a.replace('*','')+'.SDSS.astrometry.pdf' for
                 a in data2['ID']]

        out_pdf_2017a_sdss = astro_path+'GNIRS_2017A_Targets_Astrometry.SDSS.pdf'
        print '### Writing : ', out_pdf_2017a_sdss
        pdfmerge.merge(files, out_pdf_2017a_sdss)

    # + on 14/01/2017
    do_WISE = 0
    if do_WISE:
        run_main('SDF', cat2_survey='WISE')
        run_main('DEEP2', cat2_survey='WISE')

        infile2 = path0 + 'targets.2017a.txt'
        print '### Reading : ', infile2
        data2   = asc.read(infile2)
        files = [astro_path+a.replace('*','')+'.WISE.astrometry.pdf' for
                 a in data2['ID']]

        out_pdf_2017a_sdss = astro_path+'GNIRS_2017A_Targets_Astrometry.WISE.pdf'
        print '### Writing : ', out_pdf_2017a_sdss
        pdfmerge.merge(files, out_pdf_2017a_sdss)

    # + on 14/01/2017
    do_2MASS = 0
    if do_2MASS:
        run_main('SDF', cat2_survey='2MASS')
        run_main('DEEP2', cat2_survey='2MASS')

        infile2 = path0 + 'targets.2017a.txt'
        print '### Reading : ', infile2
        data2   = asc.read(infile2)
        files = [astro_path+a.replace('*','')+'.2MASS.astrometry.pdf' for
                 a in data2['ID']]

        out_pdf_2017a_sdss = astro_path+'GNIRS_2017A_Targets_Astrometry.2MASS.pdf'
        print '### Writing : ', out_pdf_2017a_sdss
        pdfmerge.merge(files, out_pdf_2017a_sdss)

#enddef

# This code has been replaced with run_main()
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

    cat_path0  = '/Users/cly/data/SDF/NBcat/SEx/catalogs/'

    infile0 = path0 + 'targets_sdf_astro.2017a.txt'
    if silent == False: print '### Reading : ', infile0
    data0   = asc.read(infile0, format='commented_header')
    print data0
    
    n_sources = len(data0)

    ID  = data0['ID']
    RA  = data0['RA']
    DEC = data0['DEC']
    c0  = coords.SkyCoord(ra=RA, dec=DEC, unit=(u.hour, u.deg))

    # Mod on 12/01/2017
    cat_files = [cat_path0+'sdf_pub2_'+f0+'.cat.fits.gz' for
                 f0 in data0['filt']]
    #cat_files = [cat_path0+f0+'emitters_allphot.fits' for f0 in data0['filt']]
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

# This code has been replaced with run_main()
def DEEP2(silent=False, verbose=True):
    '''
    Cross-match DEEP2 catalogs against SDSS to get astrometry around each
    [OIII]4363 targets from Ly et al. (2015), ApJ, 805, 45

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

    if silent == False: print '### Begin check_astrometry.DEEP2 | '+systime()

    cat_path0  = '/Users/cly/data/DEEP2/DR4/phot_cat/'

    infile0 = path0 + 'targets_deep2_astro.2017a.txt'
    if silent == False: print '### Reading : ', infile0
    data0   = asc.read(infile0, format='commented_header')
    print data0

    n_sources = len(data0)

    DEEP2_idx = [ii for ii in range(n_sources) if ('DEEP' in data0['ID'][ii])]
    print '## DEEP2_idx : ', len(DEEP2_idx)

    data0 = data0[DEEP2_idx]

    ID  = [str0.replace('*','') for str0 in data0['ID']]
    RA  = data0['RA']
    DEC = data0['DEC']
    c0  = coords.SkyCoord(ra=RA, dec=DEC, unit=(u.hour, u.deg))

    cat_files = [cat_path0+f0+'_ext.fits.gz' for f0 in data0['field']]
    print cat_files

    for ii in range(n_sources):
        if silent == False: print '### Reading : ', cat_files[ii]
        data1 = fits.getdata(cat_files[ii])
        RA1   = data1.RA_SDSS
        DEC1  = data1.DEC_SDSS

        t_c = c0[ii]
        RA_off  = (RA1 - t_c.ra.deg) * 3600.0 * \
                  np.cos(np.radians(t_c.dec.deg))
        DEC_off = (DEC1 - t_c.dec.deg) * 3600.0
        in_box = np.where((np.absolute(RA_off) <= box_size) &
                          (np.absolute(DEC_off) <= box_size))[0]
        box_data1 = data1[in_box]

        if silent == False: print '## in_box : ', len(in_box)

        out_pdf = astro_path + ID[ii]+'.astrometry.pdf'
        main(t_c, ID[ii], data1=box_data1, columns=['RA_SDSS','DEC_SDSS'],
             out_pdf=out_pdf, verbose=True)

    #endfor
    if silent == False: print '### End check_astrometry.DEEP2 | '+systime()
#enddef
