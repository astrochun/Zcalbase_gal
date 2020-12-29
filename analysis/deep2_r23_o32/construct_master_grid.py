"""
construct_master_grid
=====================

Construct the DEEP2 rest-frame master grid
"""

from chun_codes import systime

from astropy.io import fits
import numpy as np
import glob
from astropy import log


def main(path0, silent=False, verbose=True):
    """
    Combine 1-D DEEP2 spectra and de-redshift to rest-frame wavelength

    Parameters
    ----------
    path0:
      path to data (ie '/Users/cly/data/DEEP2/DR4/')

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 18 February 2018
    Edited by Reagen Leimbach August 2020
    """
    
    if not silent:
        log.info('### Begin main : '+systime())

    files_2D = glob.glob(path0 + 'DEEP2_2D_Field?.f3.fits')

    cat_files = glob.glob(path0 + 'f_current/DEEP2_Field?_all_line_fit.fits')

    n_fields = len(files_2D)

    crval1 = np.zeros(n_fields)
    cdelt1 = np.zeros(n_fields)
    naxis1 = np.zeros(n_fields)

    zspec = np.array([], dtype=np.float32)
    cdelt0 = np.array([], dtype=np.float32)
    x0_min = np.zeros(n_fields)
    x0_max = np.zeros(n_fields)

    for ii in range(n_fields):
        hdr = fits.getheader(files_2D[ii])

        crval1[ii] = hdr['CRVAL1']
        cdelt1[ii] = hdr['CDELT1']
        naxis1[ii] = hdr['NAXIS1']

        cdata = fits.getdata(cat_files[ii])
        zspec = np.append(zspec, cdata.ZSPEC)

        cdelt0 = np.append(cdelt0, cdelt1[ii]/(1+zspec))
        x0_min[ii] = np.min(cdata.LMIN0)
        x0_max[ii] = np.max(cdata.LMAX0)

    avg_cdelt = np.average(cdelt0)
    log.info('## avg_cdelt : %.3f' % avg_cdelt)

    n_pixel = np.ceil((max(x0_max)-np.min(x0_min))/avg_cdelt)
    log.info('## n_pixel : %i' % n_pixel)

    crval0 = np.min(x0_min)
    log.info('## crval0 : %.3f' % crval0)

    x0 = crval0 + avg_cdelt * np.arange(n_pixel)
    print np.min(x0), np.max(x0)
           
    if not silent:
        log.info('### End main : '+systime())
# enddef
