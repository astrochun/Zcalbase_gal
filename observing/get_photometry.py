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
data0  = asc.read(infile)
N0     = len(data0)
ID0    = data0['ID']

def SDF(silent=False, verbose=True):

    '''
    Get photometric data for SDF [OIII]4363 targets

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
    '''
    
    if silent == False: print '### Begin get_photometry.SDF | '+systime()

    SDF_path0 = '/Users/cly/data/SDF/Metallicity/Optical/2008_2014/'

    # Generator function to get just Keck and MMT targets
    SDF_idx = [ii for ii in range(N0) if ('Keck' in ID0[ii] or 'MMT' in ID0[ii])]
    data = data0[SDF_idx]
    ID   = [str0.replace('*','') for str0 in ID0[SDF_idx]]

    source_file = SDF_path0 + 'SED/source_info.sort.NB.txt'
    if silent : print '### Reading : ', source_file
    source_data = asc.read(source_file)
    source_ID   = [str0.replace('#','') for str0 in source_data['col1']]

    print ID, source_ID
    idx1, idx2 = match_nosort_str(ID, source_ID)
    print idx1, idx2
    
    SED_file = SDF_path0 + 'SED/OIII4363_sample.phot.NB.fits'
    if silent == False: print '### Reading : ', SED_file
    SED_hdu   = fits.open(SED_file)
    SED_data0 = SED_hdu[1].data

    s_idx1, s_idx2 = match_nosort(source_data['col2'][idx2], SED_data0.NUMBER)
    print s_idx1, s_idx2
    SED_data = SED_data0[s_idx2]

    filt0   = ['U','B','V','R','I','Z','J','H','K']
    s_filt0 = ["'"+f0+"'" for f0 in filt0]
    arr0    = ['mag_'+f0 for f0 in filt0]
    
    for ff in range(len(filt0)):
        cmd = 'mag_'+filt0[ff]+' = SED_data.'+filt0[ff]+'_MAG_AUTO'
        exec(cmd)

    cmd0 = "tab0 = Table([ID,"+",".join(arr0)+"],names=('ID',"+','.join(s_filt0)+"))"
    exec(cmd0) #print cmd0
    print tab0
    if silent == False: print '### End get_photometry.SDF | '+systime()
#enddef

