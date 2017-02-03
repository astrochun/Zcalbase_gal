"""
get_photometry
==============

Get photometric data for [OIII]4363 targets for Gemini-N/GNIRS 2017A program
"""

import sys, os

from chun_codes import systime
from chun_codes import match_nosort, match_nosort_str

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt
import glob

from astropy.table import Table

path0  = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
infile = path0 + 'targets.2017a.txt'
data0  = asc.read(infile, format='commented_header') # Mod on 15/01/2017
N0     = len(data0)
ID0    = data0['ID']

def SDF(mag_auto=False, silent=False, verbose=True):
    '''
    Get photometric data for SDF [OIII]4363 targets,
    Ly et al. (2016), ApJS, 226, 5

    Parameters
    ----------
    mag_auto : boolean; Default: False
      Set to use SExtractors MAG_AUTO. Otherwise, use aperture photometry with
      aperture correction (what is published in Ly et al. 2016). This is the
      FAST photometric catalog

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 11 January 2017
     - Re-organized to handle both aperture fluxes and Kron aperture fluxes
    Modified by Chun Ly, 12 January 2017
     - Re-order filter list
    Modified by Chun Ly, 2 February 2017
     - Return table of photometry
    '''
    
    if silent == False: print '### Begin get_photometry.SDF | '+systime()

    SDF_path0 = '/Users/cly/data/SDF/Metallicity/Optical/2008_2014/'

    # Generator function to get just Keck and MMT targets
    SDF_idx = [ii for ii in range(N0) if ('Keck' in ID0[ii] or
                                          'MMT' in ID0[ii])]
    data = data0[SDF_idx]
    ID   = [str0.replace('*','') for str0 in ID0[SDF_idx]]

    # Change file to the non-sorted one on 11/01/2017
    source_file = SDF_path0 + 'SED/source_info.NB.txt'
    if silent == False: print '### Reading : ', source_file
    source_data = asc.read(source_file)
    source_ID   = [str0.replace('#','') for str0 in source_data['col1']]

    idx1, idx2 = match_nosort_str(ID, source_ID)
    print '## idx1 : ', idx1
    print '## idx2 : ', idx2

    # Moved up on 11/01/2017
    filt0   = ['U','B','V','R','I','Z','J','H','K']
    s_filt0 = ["'"+f0+"'" for f0 in filt0]
    arr0    = ['mag_'+f0 for f0 in filt0]

    # Mod on 11/01/2017
    if mag_auto == True:
        SED_file = SDF_path0 + 'SED/OIII4363_sample.phot.NB.fits'
    else:
        SED_file = SDF_path0 + 'SED/OIII4363_sample.fast.NB.fits'

    if silent == False: print '### Reading : ', SED_file
    SED_hdu   = fits.open(SED_file)
    SED_data0 = SED_hdu[1].data

    s_idx1, s_idx2 = match_nosort(source_data['col2'][idx2], SED_data0.NUMBER)
    print '## s_idx1 : ', s_idx1
    print '## s_idx2 : ', s_idx2

    SED_data = SED_data0[s_idx2]
    
    for ff in range(len(filt0)):
        if mag_auto == True:
            cmd = 'mag_'+filt0[ff]+' = SED_data.'+filt0[ff]+'_MAG_AUTO'
        else:
            #cmd = 'mag_'+filt0[ff]+' = SED_data.'+filt0[ff]+'_FLUX'
            cmd = 'mag_'+filt0[ff]+' = -2.5*np.log10(SED_data.'+filt0[ff]+'_FLUX)+23.90'
        exec(cmd)

    cmd0 = "tab0 = Table([ID,"+",".join(arr0)+"], "+\
           "names=('ID',"+','.join(s_filt0)+"))"
    print cmd0
    exec(cmd0)

    print tab0
    if silent == False: print '### End get_photometry.SDF | '+systime()

    return tab0 # + on 02/02/2017
#enddef

def DEEP2(silent=False, verbose=True):

    '''
    Get photometric data for DEEP2 [OIII]4363 targets,
    Ly et al. (2015), ApJ, 805, 45

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
    Created by Chun Ly, 11 January 2017
    Modified by Chun Ly, 2 February 2017
     - Return table of photometry
    '''

    if silent == False: print '### Begin get_photometry.DEEP2 | '+systime()

    DEEP2_path0 = '/Users/cly/data/DEEP2/DR4/'

    # Generator function to get just Keck and MMT targets
    DEEP2_idx = [ii for ii in range(N0) if ('DEEP' in ID0[ii])]

    data = data0[DEEP2_idx]
    ID   = [str0.replace('*','') for str0 in ID0[DEEP2_idx]]

    # Change file to the non-sorted one on 11/01/2017
    source_file = DEEP2_path0 + 'EAZY/source_info.txt'
    if silent == False: print '### Reading : ', source_file
    source_data = asc.read(source_file)
    source_ID   = [str0.replace('2_#','') for str0 in source_data['col1']]

    idx1, idx2 = match_nosort_str(ID, source_ID)
    print '## idx1 : ', idx1
    print '## idx2 : ', idx2

    # Mod on 12/01/2017 to have a consistent order as Gemini OT
    #filt0   = ['B','R','I','SDSS_U','SDSS_G','SDSS_R','SDSS_I','SDSS_Z',
    filt0   = ['SDSS_U','B','SDSS_G','SDSS_R','R','SDSS_I','I','SDSS_Z',
               'CFHT_U','CFHT_G','CFHT_R','CFHT_I','CFHT_Z']
    s_filt0 = ["'"+f0+"'" for f0 in filt0]
    arr0    = ['mag_'+f0 for f0 in filt0]

    SED_file = DEEP2_path0 + 'EAZY/DEEP2_OIII4363_phot.fits'
    if silent == False: print '### Reading : ', SED_file
    SED_hdu   = fits.open(SED_file)
    SED_data0 = SED_hdu[1].data

    s_idx1, s_idx2 = match_nosort(source_data['col2'][idx2], SED_data0.ID)
    print '## s_idx1 : ', s_idx1
    print '## s_idx2 : ', s_idx2

    SED_data = SED_data0[s_idx2]

    for ff in range(len(filt0)):
        cmd = 'mag_'+filt0[ff]+' = -2.5*np.log10(SED_data.'+filt0[ff]+')+23.90'
        exec(cmd)

    cmd0 = "tab0 = Table([ID,"+",".join(arr0)+"], "+\
           "names=('ID',"+','.join(s_filt0)+"))"
    print cmd0
    exec(cmd0)

    print tab0
    if silent == False: print '### End get_photometry.DEEP2 | '+systime()

    return tab0 # + on 02/02/2017
#enddef

def run_all():
    '''
    Run both SDF() and DEEP2() functions from above
    '''

    # Mod on 02/02/2017
    SDF_tab0   = SDF(mag_auto=False)
    DEEP2_tab0 = DEEP2()

