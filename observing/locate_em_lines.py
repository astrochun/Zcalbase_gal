"""
locate_em_lines
===============

Generate plots illustrating location of nebular emission lines given redshift,
the OH night skyline spectrum, and the atmospheric transmission
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

from math import pi
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust

def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

def main(in_cat, out_pdf, lambda0_min, lambda0_max, R_spec,
         format='commented_header', silent=False, verbose=True):

    '''
    Main function to read in catalog of sources and output PDF file
    illustrating emission line location, OH night skyline spectrum,
    and atmospheric transmission

    Parameters
    ----------
    in_cat : string
      Filename for ASCII catalog containing source name and redshift

    out_pdf : string
      Filename to output PDF. Full path should be provided

    lambda0_min : float or double
      Minimum rest-frame wavelength for plots

    lambda0_max : float or double
      Maximum rest-frame wavelength for plots

    R_spec : float or double
      Resolution (FWHM) of spectra to convolve OH skylines with

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 17 December 2016
    '''

    co_filename = __file__

    if silent == False: print '### Begin locate_em_lines.main() | '+systime()

    if silent == False: print '### Reading : ', in_cat
    cat_data = asc.read(in_cat)

    ID    = cat_data['ID']
    zspec = cat_data['redshift']

    # Read in file with wavelengths
    # Values are in Angstroms
    lambda0_file = os.path.dirname(co_filename) + '/' + 'lambda0.txt'
    lambda0_data = asc.read(lambda0_file, format='commented_header')
    lambda0    = lambda0_data['lambda0']
    str_lines0 = lambda0_data['str_lines0']
    #[3727, 4861.32, 4958.91, 5006.84, 6562.80, 6548.10, 6583.60, 6717.42, 6730.78]
    
    rousselot_file = os.path.dirname(co_filename) + '/' + 'rousselot2000.dat'
    rousselot_data = asc.read(rousselot_file, format='commented_header')
    n_OH    = len(rousselot_data)
    n_pix   = 27000.0
    OH_bgd  = np.zeros(n_pix)
    lam_OH  = np.arange(n_pix)
    for rr in range(n_OH):
        t_lambda = rousselot_data['lambda'][rr]
        t_FWHM   = t_lambda / R_spec
        temp     = gaussian(lam_OH, t_lambda, t_FWHM/(2 * np.sqrt(2*np.log(2))))
        OH_bgd += rousselot_data['flux'][rr] * temp
    #endfor

    l_temp1 = (1+zspec) * lambda0_min
    l_temp2 = (1+zspec) * lambda0_max
    t_mark = np.where((lam_OH >= np.min(l_temp1)) & (lam_OH <= np.max(l_temp2)))[0]
    max0 = np.max(OH_bgd[t_mark])
    print '## max0 : ', max0
    
    #OII_loc   = (1.0+zspec) * 3727.00
    #Hb_loc    = (1.0+zspec) * 4861.32
    #OIIIa_loc = (1.0+zspec) * 4958.91
    #OIIIb_loc = (1.0+zspec) * 5006.84
    #Ha_loc    = (1.0+zspec) * 6562.80
    #NIIa_loc  = (1.0+zspec) * 6548.10
    #NIIb_loc  = (1.0+zspec) * 6583.60
    #SIIa_loc  = (1.0+zspec) * 6717.42
    #SIIb_loc  = (1.0+zspec) * 6730.78

    n_sources = len(ID)
    
    n_panels = 4 # Number of figure panels

    pp = PdfPages(out_pdf)

    for ii in range(n_sources):
        c_ii = ii % n_panels + 1
        if c_ii == 1: fig, (ax1, ax2, ax3, ax4) = plt.subplots(4) #, sharex=True)

        cmd1 = 'ax0 = ax'+str(c_ii)
        exec(cmd1)

        xlim = [(1.0 + zspec[ii]) * lambda0_min, (1.0 + zspec[ii]) * lambda0_max]
        ax0.set_xlim(xlim)
        ax0.set_ylim([0,1.1])
        ax0.minorticks_on()
        
        t_x0 = xlim[0] + 0.01 * (xlim[1]-xlim[0])
        ax0.annotate(ID[ii]+r' $z_{\rm spec}$ = '+str(zspec[ii]), [t_x0, 1.05],
                     xycoords='data', va='top', ha='left')

        # Overlay OH night skyline with it normalized
        ax0.plot(lam_OH, OH_bgd / max0)

        for ll in range(len(lambda0)):
            t_x = (1+zspec[ii]) * lambda0[ll]
            ax0.plot(np.repeat(t_x,2), [0,1], 'r--')
            if c_ii == 1:
                ax0.annotate(str_lines0[ll], xy=(t_x, 1.01), xycoords='data',
                             xytext=(0, 0.0), textcoords='offset points',
                             color='r', size='small', ha='center', va='bottom',
                             rotation=90)
            
        if c_ii == 4:
            ax0.set_xlabel('Wavelength (Angstroms)')
            fig.set_size_inches(8,8)
            #fig.tight_layout()
            subplots_adjust(left=0.05, bottom=0.075, top=0.975, right=0.95,
                            wspace=0.05, hspace=0.13)

            fig.savefig(pp, format='pdf') #, bbox_inches='tight')
            fig.clear()
    #endfor

    if silent == False: print '### Writing : ', out_pdf
    pp.close()

    if silent == False: print '### End locate_em_lines.main() | '+systime()
#enddef

def zcalbase_gal_gemini():

    '''
    Function to run main() but for Gemini-N/GNIRS targets

    Parameters
    ----------
    None
	  
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 17 December 2016
    '''

    print '### Begin locate_em_lines.zcalbase_gal_gemini() | '+systime()

    path0 = '/Users/cly/Dropbox/Observing/2017A/Gemini/'
    in_cat  = path0 + 'targets.txt'
    out_pdf = path0 + 'locate_em_lines.pdf'

    main(in_cat, out_pdf, 6400, 6800, 3000)
    
    print '### End locate_em_lines.zcalbase_gal_gemini() | '+systime()
#enddef
