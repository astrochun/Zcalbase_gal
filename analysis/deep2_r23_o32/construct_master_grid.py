from astropy.io import fits
import numpy as np
import glob
from os.path import join

from .log_commons import log_stdout


def main(path0, log=None):
    """
    Combine 1-D DEEP2 spectra and de-redshift to rest-frame wavelength
    Construct the DEEP2 rest-frame master grid

    :param path0: str. path to data (ie '/Users/cly/data/DEEP2/DR4/')
    :param log: LogClass or logging object

    No returns

    Notes
    -----
    Created by Chun Ly, 18 February 2018
    Edited by Reagen Leimbach August 2020
    """

    if log is None:
        log = log_stdout()

    log.info("starting ...")

    twod_files = glob.glob(join(path0, 'DEEP2_2D_Field?.f3.fits'))

    cat_files = glob.glob(join(path0,
                               'f_current/DEEP2_Field?_all_line_fit.fits'))

    n_fields = len(twod_files)

    crval1 = np.zeros(n_fields)
    cdelt1 = np.zeros(n_fields)
    naxis1 = np.zeros(n_fields)

    zspec = np.array([], dtype=np.float32)
    cdelt0 = np.array([], dtype=np.float32)
    x0_min = np.zeros(n_fields)
    x0_max = np.zeros(n_fields)

    for ii in range(n_fields):
        log.info(f"Reading: {twod_files[ii]}")
        hdr = fits.getheader(twod_files[ii])

        crval1[ii] = hdr['CRVAL1']
        cdelt1[ii] = hdr['CDELT1']
        naxis1[ii] = hdr['NAXIS1']

        log.info(f"Reading: {cat_files[ii]}")
        cat_tab = fits.getdata(cat_files[ii])
        zspec = np.append(zspec, cat_tab.ZSPEC)

        cdelt0 = np.append(cdelt0, cdelt1[ii]/(1+zspec))
        x0_min[ii] = np.min(cat_tab.LMIN0)
        x0_max[ii] = np.max(cat_tab.LMAX0)

    avg_cdelt = np.average(cdelt0)
    log.info(f"## avg_cdelt : {avg_cdelt:%.3f}")

    n_pixel = np.ceil((max(x0_max)-np.min(x0_min))/avg_cdelt)
    log.info(f"## n_pixel : {n_pixel:%i}")

    crval0 = np.min(x0_min)
    log.info(f"## crval0 : {crval0:%.3f})")

    x0 = crval0 + avg_cdelt * np.arange(n_pixel)
    log.info(np.min(x0), np.max(x0))

    log.info("finished.")
