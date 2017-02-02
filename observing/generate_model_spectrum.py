"""
generate_model_spectrum
=======================

This code generates a model spectrum based on emission-line fluxes and
expected ratios. 
"""

import sys, os

from chun_codes import systime, match_nosort_str

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import astropy.units as u

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt

from astropy.table import Table

from locate_em_lines import gaussian, gaussian_R

def main(tab0, R_spec=3000.0, out_path='', unit0=u.micron, silent=False,
         verbose=True):

    '''
    Main function() that generates a spectrum containing H-alpha, [NII],
    and [SII]

    Parameters
    ----------
    tab0 : astropy.table.table.Table
      Astropy table containing source ID, redshift, H-alpha flux, and NII/Ha
      flux ratio

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 1 February 2017
    Modified by Chun Ly, 2 February 2017
     - Minor modification for a small bug
     - Add out_path keyword option
     - Add unit0 keyword to specify desired units for wavelength
     - Fix wavelength grid to be a set wavelength in observed space
     - Minor bug fix with spectrum normalize (need dl factor)
    '''
    
    if silent == False:
        print '### Begin generate_model_spectrum.main | '+systime()

    x_min0, x_max0 = 5000.0, 7000.0 # in Angstroms
    
    n_sources = len(tab0)

    s_lambda0 = ['Ha', 'NII', 'SII', 'SII']
    lambda0 = np.array([6562.8, 6583.6, 6716.42, 6730.78])

    ID0      = tab0['ID']
    zspec    = tab0['zspec']
    Ha_flux  = 10**(tab0['Ha_flux'])
    logNIIHa = tab0['logNIIHa']
    
    scale0 = [1.0, 0.0, 0.1, 0.1]

    wave_scale = (1 * u.Angstrom).to(unit0).value

    for ii in range(n_sources):
        z0 = (1+zspec[ii])

        # Mod on 02/02/2017
        x_min, x_max = x_min0*z0, x_max0*z0
        dl = 0.25
        wave = np.arange(x_min, x_max, dl) # rest wavelength

        flux = np.zeros(len(wave))

        for ll in range(len(lambda0)):
            o_lambda = lambda0[ll] * z0
            
            s_temp = gaussian_R(wave, o_lambda, R_spec)
            if scale0[ll] != 0.0:
                flux += s_temp * Ha_flux[ii] * dl * scale0[ll]
            else:
                flux += s_temp * Ha_flux[ii] * dl * 10**(logNIIHa[ii])

        f_table = Table([wave * wave_scale, flux])

        outfile = out_path + ID0[ii]+'.spec.txt'
        if silent == False: print '### Writing : ', outfile
        asc.write(f_table, outfile, format='no_header', overwrite=True)

    if silent == False:
        print '### End generate_model_spectrum.main | '+systime()
#enddef

def gnirs_2017a(silent=False, verbose=True):
    '''
    Function to run generate_model_spectrum.main() for Gemini-N/GNIRS
    targets in 2017A program

    Parameters
    ----------
    None
          
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 2 February 2017
    '''
    
    path0   = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    infile1 = path0 + 'targets.2017a.txt'

    path1   = r'/Users/cly/Google Drive/Documents/Proposals/2017A/Gemini/'
    infile2 = path1 + 'PIT_SDF_DEEP2.bright.txt'

    if silent == False: print '### Reading : ', infile1
    tab1 = asc.read(infile1, format='commented_header')

    ID    = [str0.replace('*','') for str0 in tab1['ID']]
    zspec = tab1['redshift']

    if silent == False: print '### Reading : ', infile2
    tab2 = asc.read(infile2)

    ID2  = [str0.replace('#','') for str0 in tab2['ID']]
    idx1, idx2 = match_nosort_str(ID, ID2)

    Ha_flux  = tab2['Ha_flux'][idx2]
    logNIIHa = tab2['logNIIHa'][idx2]

    tab0 = Table([ID, zspec, Ha_flux, logNIIHa],
                 names=('ID','zspec','Ha_flux','logNIIHa'))

    print tab0

    qR_spec = 7200/2.5

    out_path = path0 + 'ETC/'
    main(tab0, R_spec, out_path=out_path, unit0=u.nm)
