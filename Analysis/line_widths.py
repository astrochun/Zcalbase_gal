"""
line_widths
====

Plot line widths from emission-line fits for DEEP2 sample
"""

from chun_codes import systime

from astropy.io import ascii as asc
from astropy.io import fits

from os.path import exists

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

    out_path = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'
    combine_stack_file = out_path + 'DEEP2_all_line_fit.fits'

    if not exists(combine_stack_file):
        for ii in range(len(files)):
            data, hdr = fits.getdata(files[ii], header=True)
            if ii == 0:
                data0 = Table(data)
                hdr0 = hdr
            else:
                data0 = vstack([data0, Table(data)])

        log.info("Writing : "+combine_stack_file)
        data0.write(combine_stack_file, format='fits')
    else:
        log.info("Reading : "+combine_stack_file)
        tab0, hdr0 = fits.getdata(combine_stack_file, header=True)
        data0 = Table(tab0)

    OIII_5007 = data0['OIIIR_FLUX_MOD'].data
    OIII_4959 = data0['OIIIB_FLUX_MOD'].data
    OII       = data0['OII_FLUX_MOD'].data
    HB        = data0['HB_FLUX_MOD'].data

    lR23 = np.log10((OII + OIII_5007+OIII_4959)/HB)
    lO32 = np.log10((OIII_5007+OIII_4959)/OII)

    out_pdf = out_path + 'line_width.pdf'
    pp = PdfPages(out_pdf)

    xlim = [0.0,10.5]

    str_lines = ['OIIB','OIIR', 'HB', 'OIIIB', 'OIIIR']
    for line in str_lines:
        fig, ax = plt.subplots(nrows=2, ncols=2)

        x_temp = data0[line+'_SIGMA'].data
        good = np.where((x_temp < 90) & (x_temp > -90))[0]
        ax[0,0].hist(x_temp[good], bins=30, alpha=0.5)

        ax[0,0].annotate(line, [0.95,0.95], xycoords='axes fraction', ha='right',
                       va='top')
        ax[0,0].set_ylabel('N', fontsize=16)
        ax[0,0].set_xticklabels([])

        ax[0,1].set_xticklabels([])

        ax[1,0].scatter(x_temp[good], lR23[good], alpha=0.5, edgecolor='none')
        ax[1,0].set_ylim(0.0,2.25)
        ax[1,0].set_xlabel(r'$\sigma$ [$\AA$]', fontsize=16)
        ax[1,0].set_ylabel(r'$\log(R_{23})$', fontsize=16)

        ax[1,1].scatter(x_temp[good], lO32[good], alpha=0.5, edgecolor='none')
        ax[1,1].set_ylim(-1,2.25)
        ax[1,1].set_xlabel(r'$\sigma$ [$\AA$]', fontsize=16)
        ax[1,1].set_ylabel(r'$\log(O_{32})$', fontsize=16)

        for t_ax in ax.flatten():
            t_ax.set_xlim(xlim)

        fig.set_size_inches(11,8.5)
        plt.subplots_adjust(left=0.06, right=0.99, top=0.97, bottom=0.065,
                            wspace=0.17, hspace=0.01)
        fig.savefig(pp, format='pdf')
    #endfor

    log.info('Writing : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

