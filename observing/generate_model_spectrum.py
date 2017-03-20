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
import astropy.constants as const

import numpy as np
import array, time, sets

import matplotlib.pyplot as plt

from astropy.table import Table

from locate_em_lines import gaussian, gaussian_R
import get_photometry

def main(tab0, velocity=100.0*u.km/u.s, out_path='', unit0=u.micron,
         silent=False, verbose=True):

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
     - Changed R_spec to km/s to provide an intrinsic gas velocity
       for nebular emission lines. Most ETC will factor in spectral
       resolution
     - Adopt non-zero continuum if magnitude column ('mag') is present in table
    Modified by Chun Ly, 2 February 2017
     - Remove lambda_cen keyword option
     - Assume continuum is flat in F-nu
    '''
    
    if silent == False:
        print '### Begin generate_model_spectrum.main | '+systime()

    x_min0, x_max0 = 5000.0, 9000.0 # in Angstroms
    
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
        dl = 0.50 # in Angstrom
        wave = np.arange(x_min, x_max, dl) # rest wavelength

        flux = np.zeros(len(wave))

        # Normalize spectrum based on optical data | + on 02/02/2017
        if 'mag' in tab0.colnames:
            f_nu = 10**(-0.4*(tab0['mag'][ii]-23.90)) * u.microJansky
            f_nu = f_nu.to(u.erg/u.s/(u.cm**2)/u.Hz)
            #Assume that spectrum is flat in F_nu
            f_lambda = f_nu * const.c / ((wave*u.Angstrom)**2)
            f_lambda0 = f_lambda.to(u.erg/u.s/(u.cm**2)/u.Angstrom).value
            #print ii, f_nu, f_lambda0
            flux = f_lambda0 #np.repeat(f_lambda0, len(wave))

        for ll in range(len(lambda0)):
            o_lambda = lambda0[ll] * z0

            # This is actually the inverse. Equivalent to a spectral resolution
            # since gaussian_R accepts R_spec
            R_vel = 1.0/(velocity/const.c.to(u.km/u.s))
            #print ID0[ii], s_lambda0[ll], R_vel
            s_temp = gaussian_R(wave, o_lambda, R_vel)
            s_temp /= np.sum(s_temp)
            #if ll == 0: print ii, np.sum(s_temp)
            if scale0[ll] != 0.0:
                flux += s_temp * Ha_flux[ii] * scale0[ll] / dl
            else:
                flux += s_temp * Ha_flux[ii] * 10**(logNIIHa[ii]) / dl

        f_table = Table([wave * wave_scale, flux])

        outfile = out_path + ID0[ii]+'.spec.sed' # Mod on 20/03/2017
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

    if silent == False: print '### Reading : ', infile1
    tab1 = asc.read(infile1, format='commented_header')

    ID    = [str0.replace('*','') for str0 in tab1['ID']]
    zspec = tab1['redshift']

    path1   = r'/Users/cly/Google Drive/Documents/Proposals/2017A/Gemini/'
    infile2 = path1 + 'PIT_SDF_DEEP2.bright.txt'

    if silent == False: print '### Reading : ', infile2
    tab2 = asc.read(infile2)

    ID2  = [str0.replace('#','') for str0 in tab2['ID']]
    idx1, idx2 = match_nosort_str(ID, ID2)

    Ha_flux  = tab2['Ha_flux'][idx2]
    logNIIHa = tab2['logNIIHa'][idx2]

    # Get photometry
    # later + on 02/02/2017
    SDF_phot   = get_photometry.SDF(verbose=False)
    DEEP2_phot = get_photometry.DEEP2(verbose=False)

    s_idx1, s_idx2 = match_nosort_str(ID, SDF_phot['ID'])
    d_idx1, d_idx2 = match_nosort_str(ID, DEEP2_phot['ID'])

    #t0 = Table([DEEP2_phot['SDSS_Z'], DEEP2_phot['CFHT_Z'], DEEP2_phot['I']])
    #print t0

    i_mag1 = np.where(np.isfinite(DEEP2_phot['SDSS_Z']))[0]
    i_mag2 = np.where(np.isfinite(DEEP2_phot['CFHT_Z']))[0]
    deep2_mag = DEEP2_phot['I']
    deep2_mag[i_mag1] = DEEP2_phot['SDSS_Z'][i_mag1]
    deep2_mag[i_mag2] = DEEP2_phot['CFHT_Z'][i_mag2]

    mag_phot0         = np.zeros(len(ID))
    mag_phot0[s_idx1] = SDF_phot['J'][s_idx2]

    mag_phot0[d_idx1] = deep2_mag[d_idx2]

    print mag_phot0

    tab0 = Table([ID, zspec, Ha_flux, logNIIHa, mag_phot0],
                 names=('ID','zspec','Ha_flux','logNIIHa','mag'))

    print tab0

    vel0 = 100.0 * u.km/u.s

    out_path = path0 + 'ETC/'
    main(tab0, velocity=vel0, out_path=out_path, unit0=u.nm)
    #lambda_cen=1.25*u.micron)

def gnirs_2017b(silent=False, verbose=True):
    '''
    Function to run generate_model_spectrum.main() for Gemini-N/GNIRS
    targets in 2017B program (pre-proposal planning)

    Parameters
    ----------
    None

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 20 March 2017
    '''

    path0   = r'/Users/cly/Google Drive/Documents/Proposals/2017B/Gemini/'
    infile1 = path0 + 'targets.2017b.txt'

    if silent == False: print '### Reading : ', infile1
    tab1 = asc.read(infile1, format='commented_header')

    ID    = [str0.replace('*','') for str0 in tab1['ID']]
    zspec = tab1['redshift']

    path1   = r'/Users/cly/Google Drive/Documents/Proposals/2017A/Gemini/'
    infile2 = path1 + 'PIT_SDF_DEEP2.bright.txt'

    if silent == False: print '### Reading : ', infile2
    tab2 = asc.read(infile2)

    ID2  = [str0.replace('#','') for str0 in tab2['ID']]
    idx1, idx2 = match_nosort_str(ID, ID2)

    Ha_flux  = tab2['Ha_flux'][idx2]
    logNIIHa = tab2['logNIIHa'][idx2]

    # Get photometry
    # later + on 02/02/2017
    SDF_phot   = get_photometry.SDF(verbose=False)
    DEEP2_phot = get_photometry.DEEP2(verbose=False)

    print DEEP2_phot

    s_idx1, s_idx2 = match_nosort_str(ID, SDF_phot['ID'])
    d_idx1, d_idx2 = match_nosort_str(ID, DEEP2_phot['ID'])

    #t0 = Table([DEEP2_phot['SDSS_Z'], DEEP2_phot['CFHT_Z'], DEEP2_phot['I']])
    #print t0

    i_mag1 = np.where(np.isfinite(DEEP2_phot['SDSS_Z']))[0]
    i_mag2 = np.where(np.isfinite(DEEP2_phot['CFHT_Z']))[0]
    deep2_mag = DEEP2_phot['I']
    deep2_mag[i_mag1] = DEEP2_phot['SDSS_Z'][i_mag1]
    deep2_mag[i_mag2] = DEEP2_phot['CFHT_Z'][i_mag2]

    mag_phot0         = np.zeros(len(ID))
    mag_phot0[s_idx1] = SDF_phot['J'][s_idx2]

    mag_phot0[d_idx1] = deep2_mag[d_idx2]

    print mag_phot0

    tab0 = Table([ID, zspec, Ha_flux, logNIIHa, mag_phot0],
                 names=('ID','zspec','Ha_flux','logNIIHa','mag'))

    print tab0

    vel0 = 100.0 * u.km/u.s

    out_path = path0 + 'ETC/'
    main(tab0, velocity=vel0, out_path=out_path, unit0=u.nm)
    #lambda_cen=1.25*u.micron)
