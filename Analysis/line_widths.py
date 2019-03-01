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
from matplotlib.backends.backend_pdf import PdfPages

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

    OIII_5007 = data0['OIIIR_FLUX_MOD'].data
    OIII_4959 = data0['OIIIB_FLUX_MOD'].data
    OII       = data0['OII_FLUX_MOD'].data
    HB        = data0['HB_FLUX_MOD'].data

    lR23 = np.log10((OII + OIII_5007+OIII_4959)/HB)
    lO32 = np.log10((OIII_5007+OIII_4959)/OII)

    out_pdf = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'+\
              'line_width.pdf'
    pp = PdfPages(out_pdf)

    str_lines = ['OIIB','OIIR', 'HB', 'OIIIB', 'OIIIR']
    for line in str_lines:
        fig, ax = plt.subplots(nrows=2)

        x_temp = data0[line+'_SIGMA'].data
        good = np.where((x_temp < 90) & (x_temp > -90))[0]
        ax[0].hist(x_temp[good], bins=30, alpha=0.5)

        ax[0].annotate(line, [0.95,0.95], xycoords='axes fraction', ha='right',
                       va='top')

        ax[1].scatter(x_temp[good], lR23[good], alpha=0.5, edgecolor='none')
        ax[1].set_ylim(0.0,2.0)
        ax[1].set_xlabel(r'$\sigma$ [$\AA$]')
        ax[1].set_ylabel(r'$\log(R_{23})$')
        fig.savefig(pp, format='pdf')
    #endfor

    log.info('Writing : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

