"""
local_analog_calibration
====

Plots R23 and O32 line ratios and compare it to Bian et al. (2018)
calibration that is based on local analogs to z~2 galaxies
"""

import sys, os

from chun_codes import systime, match_nosort_str

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

from .green_peas_calibration import get_zcalbase_sample

def bian18_R23(OH):
    '''
    Function to return log(R23) given metallicity

    Parameters
    ----------
     OH : 12+log(O/H)
    '''

    R23_coeff = [-0.32293, 7.2954, -54.8284, 138.0430]
    R23_p = np.poly1d(R23_coeff)

    return R23_p(OH)
#enddef

def bian18_O32(O32):
    '''
    Function to return metallicity given log(O32)

    Parameters
    ----------
     O32 : log([OIII]/[OII])
    '''

    OH = 8.54 - 0.59 * O32

    return OH
#enddef

def main(lR23, lO32, OH, out_pdf, ID=[], lR23_err=[], lO32_err=[],
         OH_err=[], R23_xra=[], O32_xra=[], yra=[], ctype=[], label=[''],
         silent=False, verbose=True):

    '''
    Main function to plot dataset against Bian+ (2018) calibration

    Parameters
    ----------
    lR23 : log(R_23)
    lO32 : log([OIII]/[OII])
    OH   : 12+log(O/H)

    out_pdf : full path for output PDF

    ID   : ID of source
      Provide if plots should be annotated

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 29 November 2018
    '''

    if silent == False: log.info('### Begin main : '+systime())

    fig, ax = plt.subplots(ncols=2)

    n_sample = len(lR23)

    min1, max1 = np.zeros(n_sample), np.zeros(n_sample)
    OH_min1, OH_max1 = np.zeros(n_sample), np.zeros(n_sample)
    for nn in range(n_sample):
        min1[nn] = np.min(lO32[nn])
        max1[nn] = np.max(lO32[nn])

        OH_min1[nn] = np.min(OH[nn])
        OH_max1[nn] = np.max(OH[nn])

    O32_min = np.min(min1)
    O32_max = np.max(max1)
    O32_arr = np.arange(O32_min,O32_max, 0.025)

    bian_OH = bian18_O32(O32_arr)

    # Grid of 12+log(O/H)
    if len(yra) == 0:
        print(len(yra))
        print(min(OH_min1))
        OH_max = np.nanmax(OH_max1)
        print(OH_max)
        OH_arr = np.arange(min(OH_min1),OH_max,0.25)
    else:
        OH_arr = np.arange(yra[0],yra[1],0.05)
    print('OH_arr:', OH_arr)
    bian_R23 = bian18_R23(OH_arr)

    if len(ctype) == 0:
        ctype = ['blue'] * n_sample
    if len(ctype) == 0:
        label = [''] * n_sample

    for nn in range(n_sample):
        ax[0].scatter(lR23[nn], OH[nn], color=ctype[nn], edgecolor='none',
                      marker='o', alpha=0.5, label=label[nn])
        if len(ID) != 0:
            for ii in range(len(lR23[nn])):
                ax[0].annotate(ID[nn][ii], [lR23[nn][ii], OH[nn][ii]],
                               xycoords='data', size=4, va='bottom',
                               ha='left')

        if len(OH_err) != 0:
            ax[0].errorbar(lR23[nn], OH[nn], yerr=OH_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5,
                           fmt=None, label=None)

        if len(lR23_err) != 0:
            ax[0].errorbar(lR23[nn], OH[nn], xerr=lR23_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5,
                           fmt=None, label=None)

    ax[0].plot(bian_R23, OH_arr, 'k--', label='Bian+(2018)')

    if len(R23_xra) != 0: ax[0].set_xlim(R23_xra)
    if len(yra) != 0: ax[0].set_ylim(yra)

    ax[0].set_xlabel(r'$\log(R_{23})$')
    ax[0].set_ylabel(r'$12+\log({\rm O/H})$')

    ax[0].legend(loc='lower left', framealpha=0.5, fontsize=10)

    for nn in range(n_sample):
        ax[1].scatter(lO32[nn], OH[nn], color=ctype[nn], edgecolor='none',
                      marker='o', alpha=0.5)
        if len(ID) != 0:
            for ii in range(len(lO32[nn])):
                ax[1].annotate(ID[nn][ii], [lO32[nn][ii], OH[nn][ii]],
                               xycoords='data', size=4, va='bottom',
                               ha='left')

        if len(OH_err) != 0:
            ax[1].errorbar(lO32[nn], OH[nn], yerr=OH_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5,
                           fmt=None, label=None)

        if len(lO32_err) != 0:
            ax[1].errorbar(lO32[nn], OH[nn], xerr=lO32_err[nn], mec=ctype[nn],
                           ecolor=ctype[nn], capsize=0, alpha=0.5,
                           fmt=None, label=None)

    ax[1].plot(O32_arr, bian_OH, 'k--', label='Bian+(2018)')

    if len(O32_xra) != 0: ax[1].set_xlim(O32_xra)
    if len(yra) != 0: ax[1].set_ylim(yra)

    ax[1].set_yticklabels([])
    ax[1].set_xlabel(r'$\log(O_{32})$')

    plt.subplots_adjust(left=0.065, right=0.99, bottom=0.1, top=0.99,
                        wspace=0.01)
    fig.set_size_inches(10,5)
    fig.savefig(out_pdf,format='pdf')

    if silent == False: log.info('### End main : '+systime())
#enddef

def get_measurements(data):

    lR23 = np.log10(data['R23'].data)
    lO32 = np.log10(data['O32'].data)
    OH   = data['OH'].data

    OH_err   = np.row_stack((data['OH_lo'].data,data['OH_hi'].data))

    lR23_lo = lR23 - np.log10(data['R23'].data - data['R23_lo'].data)
    lR23_hi = np.log10(data['R23'].data + data['R23_hi'].data) - lR23
    lR23_err = np.row_stack((lR23_lo,lR23_hi))

    lO32_lo = lO32 - np.log10(data['O32'].data - data['O32_lo'].data)
    lO32_hi = np.log10(data['O32'].data + data['O32_hi'].data) - lO32
    lO32_err = np.row_stack((lO32_lo,lO32_hi))

    return lR23, lO32, OH, OH_err, lR23_err, lO32_err
#enddef

def get_DEEP2(path0):
    infile = path0 + 'DEEP2_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_measurements(data)

    return data, lR23, lO32, OH, OH_err, lR23_err, lO32_err
#enddef

def get_MACT(path0):
    infile = path0 + 'MACT_R23_O32_derived.tbl'

    log.info('### Reading : '+infile)

    data = asc.read(infile)

    lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_measurements(data)

    return data, lR23, lO32, OH, OH_err, lR23_err, lO32_err
#enddef

def DEEP2_OIII4363():
    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_DEEP2(path0)

    ID = data['ID']

    out_pdf = path0 + 'DEEP2_R23_O32_Bian18.pdf'
    main([lR23], [lO32], [OH], out_pdf, ID=[ID], lR23_err=[lR23_err],
         lO32_err=[lO32_err], OH_err=[OH_err], R23_xra=[0.75,1.05],
         O32_xra=[0.05,0.95], yra=[7.1,8.65], label=[r'DEEP2 (Ly+2015)'])

#enddef

def MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data, lR23, lO32, OH, OH_err, lR23_err, lO32_err = get_MACT(path0)

    ID = data['ID']

    out_pdf = path0 + 'MACT_R23_O32_Bian18.pdf'
    main([lR23], [lO32], [OH], out_pdf, ID=[ID], lR23_err=[lR23_err],
         lO32_err=[lO32_err], OH_err=[OH_err], R23_xra=[0.6,1.15],
         O32_xra=[-0.55,2.1], yra=[7.1,8.85],
         label=[r'$\mathcal{MACT}$  (Ly+2016)'])


#enddef

def DEEP2_MACT_OIII4363():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    # DEEP2
    DEEP2_data, DEEP2_lR23, DEEP2_lO32, DEEP2_OH, DEEP2_OH_err, \
        DEEP2_lR23_err, DEEP2_lO32_err = get_DEEP2(path0)

    # MACT
    MACT_data, MACT_lR23, MACT_lO32, MACT_OH, MACT_OH_err, \
        MACT_lR23_err, MACT_lO32_err = get_MACT(path0)

    # Combine together
    lR23 = [DEEP2_lR23, MACT_lR23]
    lO32 = [DEEP2_lO32, MACT_lO32]
    OH   = [DEEP2_OH,   MACT_OH]

    OH_err = [DEEP2_OH_err, MACT_OH_err]

    lR23_err = [DEEP2_lR23_err, MACT_lR23_err]
    lO32_err = [DEEP2_lO32_err, MACT_lO32_err]

    ID = [DEEP2_data['ID'].data, MACT_data['ID'].data]

    out_pdf = path0 + 'DEEP2_MACT_R23_O32_Bian18.pdf'

    labels = [r'DEEP2 [OIII]$\lambda$4363-detected',
              r'$\mathcal{MACT}$  (Ly+2016)']
    main(lR23, lO32, OH, out_pdf, ID=ID, lR23_err=lR23_err, lO32_err=lO32_err,
         OH_err=OH_err, R23_xra=[0.6,1.15], O32_xra=[-0.55,2.1],
         yra=[7.1,8.85], ctype=['blue','green'],
         label=labels)

    out_pdf = path0 + 'DEEP2_MACT_R23_O32_Bian18.nolabel.pdf'

    main(lR23, lO32, OH, out_pdf, lR23_err=lR23_err, lO32_err=lO32_err,
         OH_err=OH_err, R23_xra=[0.6,1.15], O32_xra=[-0.55,2.1],
         yra=[7.1,8.85], ctype=['blue','green'],
         label=labels)
#enddef

def Zcalbase():

    path0 = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    ref_name0 = ['Berg2012','Kennicutt2003','Izotov1994','Thuan1995',
                 'Izotov1997','Guseva2009', 'Izotov2012', 'SDSS']
    dir0      = ['','','BCGs','BCGs','Pilyugin2012/Izotov1997',
                 'Pilyugin2012/Guseva2009', 'Pilyugin2012/Izotov2012','']

    lR23_all = []
    lO32_all = []
    OH_all   = []

    for name,dir in zip(ref_name0,dir0):
        lR23, lO32, OH = get_zcalbase_sample(name, dir=dir)

        lR23_all.append(lR23)
        lO32_all.append(lO32)
        OH_all.append(OH)

    out_pdf = path0 + 'Zcalbase_Bian18.pdf'
    label = ref_name0 #['Kennicutt+2003']
    main(lR23_all, lO32_all, OH_all, out_pdf, R23_xra=[0.1,1.05],
         O32_xra=[-0.5,1.5], yra=[7.0,9.0], label=label,
         ctype=['blue','cyan','green','yellow','red','magenta','gray','black'])
