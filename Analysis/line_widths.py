"""
line_widths
====

Plot line widths from emission-line fits for DEEP2 sample
"""

from chun_codes import systime

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from glob import glob

from astropy.table import Table, vstack
from astropy import log

def main(silent=False, verbose=True):

    '''
    Main function to plot line widths of emission-line fits

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
    Created by Chun Ly, 25 February 2019
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    path0 = '/Users/cly/data/DEEP2/DR4/f_current/'
    files = glob(path0+'*all_line_fit.fits')

    for ii in range(len(files)):
        data, hdr = fits.getdata(files[ii], header=True)
        if ii == 0:
            data0 = Table(data)
            hdr0 = hdr
        else:
            data0 = vstack([data0, Table(data)])


    str_lines = ['OIIB','OIIR', 'HB', 'OIIIB', 'OIIIR']
    for line in str_lines:
        x_temp = data0[line+'_SIGMA'].data
        good = np.where((x_temp < 90) & (x_temp > -90))[0]
        plt.hist(x_temp[good], bins=30, alpha=0.5)

    if silent == False: log.info('### End main : '+systime())
#enddef

