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

from scipy.interpolate import interp1d # + on 17/12/2016

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust

import astropy.units as u # + on 04/01/2017

import string # + on 04/01/2017

# Mauna Kea atmospheric transmission
# + on 17/12/2016
# Moved up on 21/12/2016
co_filename = __file__
atmo_file = os.path.dirname(co_filename) + '/' + 'mktrans_zm_10_10.dat'
atmo_data = asc.read(atmo_file)

# Read in file with wavelengths
# Values are in Angstroms
# Moved up on 21/12/2016
lambda0_file = os.path.dirname(co_filename) + '/' + 'lambda0.txt'
lambda0_data = asc.read(lambda0_file, format='commented_header')
lambda0    = lambda0_data['lambda0']
str_lines0 = lambda0_data['str_lines0']

def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
#enddef

def gaussian_R(x_arr, lambda_cen, R_spec):
    '''
    Generate an array consisting of a Gaussian profile given the
    spectral resolution
    
    Parameters
    ----------
    x_arr : array
      An array of wavelengths

    x_lambda : float
      The central wavelength of the Gaussian line

    R_spec : float or double
      Spectral resolution to consider width of emission lines (e.g., R = 3000)

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 20 December 2016
    '''
    
    t_FWHM = lambda_cen / R_spec # FWHM based on the wavelength of interest
    temp   = gaussian(x_arr, lambda_cen, t_FWHM/(2 * np.sqrt(2*np.log(2))))
    return temp
#enddef

def atmo_trans_val(lambda_val, R_spec, silent=True, verbose=False):
    '''
    Function to compute atmospheric transmission at given wavelength,
    factors in spectral resolution

    Parameters
    ----------
    lambda_val : float or array
      Observed wavelength of nebular emission lines

    R_spec : float or double
      Spectral resolution to consider width of emission lines (e.g., R = 3000)

    Returns
    -------
    trans0 : array
      Transmission percentage for provided emission lines

    Notes
    -----
    Created by Chun Ly, 17 December 2016
    Modified by Chun Ly, 20 December 2016
     - Weigh transmission by location on the emission line
    '''

    f_atmo = interp1d(atmo_data['col1'] * 1E4, atmo_data['col2'])

    #trans0 = f_atmo(lambda_val)
    trans0 = np.zeros(len(lambda_val))

    # Weigh transmission by location on the emission line | + on 20/12/2016
    for ii in range(len(lambda_val)):
        x_temp = np.arange(lambda_val[ii]-10,lambda_val[ii]+1,0.1)
        scale  = gaussian_R(x_temp, lambda_val[ii], R_spec)
        y_temp = f_atmo(x_temp)
        trans0[ii] = np.sum(y_temp * scale) / np.sum(scale)
    print trans0

    trans0 *= 100

    s_trans0 = ['%5.2f' % trans0[xx] for xx in range(len(trans0))]
    return trans0, s_trans0
#enddef

def OH_contam_val(OH_bgd, OH_max0, lambda_OH, lambda_val, R_spec,
                  silent=True, verbose=False):
    '''
    Function to compute OH night skyline contamination at given wavelength,
    factors in spectral resolution

    Parameters
    ----------
    OH_bgd : astropy Table
      OH skyline data (Get from Gemini webpage) 

    OH_max0 : float
      Maximum value within a certain wavelength range to normalize the
      OH skyline data to

    lambda_OH : float or array
      Wavelength of OH skyline spectrum

    R_spec : float or double
      Spectral resolution to consider width of emission lines (e.g., R = 3000)

    Returns
    -------
    contam0 : array
      OH skyline contamination for each value in lambda_val

    Notes
    -----
    Created by Chun Ly, 21 December 2016
    '''

    f_OH = interp1d(lambda_OH, OH_bgd / OH_max0)

    contam0 = np.zeros(len(lambda_val))

    # Weigh OH skyline by location on the emission line
    for ii in range(len(lambda_val)):
        x_temp = np.arange(lambda_val[ii]-10,lambda_val[ii]+1,0.1)
        scale  = gaussian_R(x_temp, lambda_val[ii], R_spec)
        y_temp = f_OH(x_temp)
        contam0[ii] = np.sum(y_temp * scale) / np.sum(scale)

    s_contam0 = ['%5.2f' % contam0[xx] for xx in range(len(contam0))]
    return contam0, s_contam0
#enddef

def overlay_filter_trans(instrument, lambda_val, in_range, ax0=None,
                         silent=True, verbose=False):
    '''
    Overlay filter transmission curves for different telescope/instrument

    Parameters
    ----------
    instrument : string
      Name of instrument. Accepted values are currently: 'GNIRS'
      The ASCII files should be placed in Filters/TELE/[instrument].
      See __filename__ for full path

    Returns
    -------
    None.

    Notes
    -----
    Created by Chun Ly, 4 January 2017
    Modified by Chun Ly, 5 January 2017
     - Fix bug in truncation of str_annot for ax0.annotate()
    '''

    if silent == False:
        print '### Begin locate_em_lines.overlay_filter_trans() | '+systime()

    if instrument == 'GNIRS':
        trans_path = os.path.dirname(co_filename) + '/Filters/Gemini/GNIRS/'
        filts  = ['X','J','H','K']
        files  = [trans_path+a.lower()+'_bl.dat' for a in filts]
        x_unit = u.micron # Unit of the wavelength
        
    if silent == False: print '### [files] : ', files

    if ax0 == None: ax0 = plt.gca()
    xlim = ax0.get_xlim()

    f_trans_val = []
    str_annot   = [a+': \n' for a in str_lines0.data[in_range]]

    for ii in range(len(files)):
        if silent == False: print '### Reading : ', files[ii]
        data    = asc.read(files[ii])
        x_scale = x_unit.to(u.angstrom)
        x_Ang   = data['col1']*x_scale
        in_plot = np.where(((xlim[0] > np.min(x_Ang)) &
                            (xlim[0] < np.max(x_Ang))) | 
                           ((xlim[1] > np.min(x_Ang)) &
                            (xlim[1] > np.max(x_Ang))))[0]

        if len(in_plot) > 0:
            ax0.plot(x_Ang, data['col2'], '--', linewidth=1.0,
                     label=instrument+' '+filts[ii])

            # Get transmission for each wavelength
            f_trans = interp1d(x_Ang, data['col2'])
            f_trans_val.append(f_trans(lambda_val))

            str_annot = [a.replace('\n','')+('%5.2f' % b)+'\n' for
                         a,b in zip(str_annot,f_trans_val[-1])]
        #endif
    #endfor
    #str_anot = [a+'\n' for a in str_annot]

    #alpha_str = list(string.lowercase)
    #str_annot = [a+str(b) for a,b in zip(str_annot, f_trans_val[-1])]

    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", alpha=0.75,
                      ec="k", lw=0.5)
    t_x0 = xlim[1] - 0.02 * (xlim[1]-xlim[0])
    # Fix bug on 05/01/2017 for truncation
    ax0.annotate(r''.join(str_annot)[:-1], [t_x0, 0.9], fontsize='x-small',
                 xycoords='data', va='top', ha='right', alpha=0.5,
                 bbox=bbox_props)

    # Legend showing the transmission filter names
    ax0.legend(loc='lower right', bbox_to_anchor=(0.995,0.03),
               fontsize='small', framealpha=0.75)

    if silent == False:
        print '### End locate_em_lines.overlay_filter_trans() | '+systime()
#enddef

def main(in_cat, out_pdf, lambda0_min, lambda0_max, R_spec, instrument,
         format='commented_header', silent=False, verbose=True):

    '''
    Main function to read in catalog of sources and output PDF file
    illustrating emission line location, OH night skyline spectrum,
    and atmospheric transmission

    Parameters
    ----------
    in_cat : string
      Filename for ASCII catalog containing source name and redshift
      in two separate columns. These should be named 'ID' and 'redshift'

    out_pdf : string
      Filename to output PDF. Full path should be provided

    lambda0_min : float or double
      Minimum rest-frame wavelength for plots.
      Note: X-axes will be in observed wavelengths

    lambda0_max : float or double
      Maximum rest-frame wavelength for plots
      Note: X-axes will be in observed wavelengths

    R_spec : float or double
      Spectral resolution to convolve OH skylines with (e.g., R = 3000)

    format : string
      Format of in_cat ASCII file. Default: "commented_header"

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------
    
    Notes
    -----
    Created by Chun Ly, 17 December 2016
    Modified by Chun Ly, 17 December 2016
     - Added atmospheric transmission
     - Fix bug with last page not being plotted. Last page now outputted on
       last source
     - Fix plotting to avoid panels being shown when sources are not available
    Modified by Chun Ly, 17 December 2016
     - Call gaussian_R() instead of using 
    Modified by Chun Ly, 4 January 2017
     - Add call to overlay_filter_trans()
    '''

    if silent == False: print '### Begin locate_em_lines.main() | '+systime()

    if silent == False: print '### Reading : ', in_cat
    cat_data = asc.read(in_cat)

    ID    = cat_data['ID']
    zspec = cat_data['redshift']

    #OII_loc   = (1.0+zspec) * 3727.00
    #Hb_loc    = (1.0+zspec) * 4861.32
    #OIIIa_loc = (1.0+zspec) * 4958.91
    #OIIIb_loc = (1.0+zspec) * 5006.84
    #Ha_loc    = (1.0+zspec) * 6562.80
    #NIIa_loc  = (1.0+zspec) * 6548.10
    #NIIb_loc  = (1.0+zspec) * 6583.60
    #SIIa_loc  = (1.0+zspec) * 6717.42
    #SIIb_loc  = (1.0+zspec) * 6730.78

    # OH night skyline file from Rousselot et al. (2000)
    rousselot_file = os.path.dirname(co_filename) + '/' + 'rousselot2000.dat'
    rousselot_data = asc.read(rousselot_file, format='commented_header')
    n_OH    = len(rousselot_data)
    n_pix   = 27000.0
    OH_bgd  = np.zeros(n_pix)
    lambda_OH  = np.arange(n_pix)
    for rr in range(n_OH):
        t_lambda = rousselot_data['lambda'][rr]
        temp     = gaussian_R(lambda_OH, t_lambda, R_spec) # Mod on 20/12/2016
        OH_bgd  += rousselot_data['flux'][rr] * temp
    #endfor

    # Get maximum OH skyline to normalize spectrum for all panels
    l_temp1 = (1+zspec) * lambda0_min
    l_temp2 = (1+zspec) * lambda0_max
    t_mark  = np.where((lambda_OH >= np.min(l_temp1)) & (lambda_OH <= np.max(l_temp2)))[0]
    OH_max0 = np.max(OH_bgd[t_mark])
    print '## OH_max0 : ', OH_max0


    n_sources = len(ID)
    print '## n_sources : ', n_sources
    
    # Number of figure panels on a page. Do not change, will break code below
    n_panels = 4 

    pp = PdfPages(out_pdf)

    for ii in range(n_sources):
        c_ii = (ii % n_panels) + 1
        # Mod on 17/12/2016 to simplify ax_arr
        if c_ii == 1: fig, ax_arr = plt.subplots(n_panels) #, sharex=True)

        ax0 = ax_arr[c_ii-1]

        #if c_ii == 1: fig, (ax1, ax2, ax3, ax4) = plt.subplots(4) #, sharex=True)
        #
        #cmd1 = 'ax0 = ax'+str(c_ii)
        #exec(cmd1)

        xlim = [(1.0 + zspec[ii]) * lambda0_min, (1.0 + zspec[ii]) * lambda0_max]
        ax0.set_xlim(xlim)
        ax0.set_ylim([0,1.2])
        ax0.minorticks_on()
        
        t_x0 = xlim[0] + 0.01 * (xlim[1]-xlim[0])
        ax0.annotate(ID[ii]+r' $z_{\rm spec}$ = '+str(zspec[ii]), [t_x0, 1.15],
                     xycoords='data', va='top', ha='left')

        # Overlay OH night skyline with it normalized
        ax0.plot(lambda_OH, OH_bgd / OH_max0)

        # Overlay atmospheric transmission | + on 17/07/2016
        ax0.plot(atmo_data['col1'].data*1e4, atmo_data['col2'].data, 'k',
                 linewidth=0.5)

        # Moved up on 04/01/2017
        in_range   = np.where((lambda0 >= lambda0_min) &
                              (lambda0 <= lambda0_max))[0]
        lambda_val = lambda0[in_range] * (1+zspec[ii])

        # Overlay filter transmission | + on 04/01/2017
        overlay_filter_trans(instrument, lambda_val, in_range, ax0=ax0) #, silent=False)
        
        # Draw emission lines
        for ll in range(len(lambda0)):
            t_x = (1+zspec[ii]) * lambda0[ll]
            ax0.plot(np.repeat(t_x,2), [0,1], 'r--')
            if c_ii == 1:
                ax0.annotate(str_lines0[ll], xy=(t_x, 1.01), xycoords='data',
                             xytext=(0, 0.0), textcoords='offset points',
                             color='r', size='small', ha='center', va='bottom',
                             rotation=90)

        # Get atmospheric transmission at lines | + on 17/12/2016
        # Mod on 20/12/2016
        trans0, s_trans0 = atmo_trans_val(lambda_val, R_spec)

        contam0, s_contam0 = OH_contam_val(OH_bgd, OH_max0, lambda_OH,
                                           lambda_val, R_spec)

        # Annotated plot | + on 20/12/2016
        str_annot = [a+': '+b+', '+c+'\n' for
                     a,b,c in zip(str_lines0.data[in_range], s_trans0,
                                  s_contam0)]

        bbox_props = dict(boxstyle="square,pad=0.3", fc="white", alpha=0.75,
                          ec="k", lw=0.5)
        
        t_x0 = xlim[0] + 0.02 * (xlim[1]-xlim[0])
        ax0.annotate(r''.join(str_annot)[:-2], [t_x0, 0.9], fontsize='x-small',
                     xycoords='data', va='top', ha='left', alpha=0.5,
                     bbox=bbox_props)
        
        # + on 17/12/2016
        if (ii == n_sources-1) and (n_sources % n_panels != 0):
            for aa in np.arange(n_sources % n_panels, n_panels):
                ax_arr[aa].axis('off')

        # Bug found on 17/12/2016. Last page excluded. Require output on last source
        if c_ii == n_panels or ii == n_sources-1: 
            ax0.set_xlabel('Wavelength (Angstroms)')
            fig.set_size_inches(8,8)
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

    R_spec = 3000.0 # Resolution of spectrograph
    # main(in_cat, out_pdf, 6400, 6800, R_spec, 'GNIRS')

    # + on 04/01/2017
    in_cat2  = path0 + 'targets.2017a.txt'
    out_pdf2 = path0 + 'locate_em_lines.2017a.pdf'
    main(in_cat2, out_pdf2, 6400, 6800, R_spec, 'GNIRS')

    print '### End locate_em_lines.zcalbase_gal_gemini() | '+systime()
#enddef
