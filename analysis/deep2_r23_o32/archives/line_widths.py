"""
line_widths
====

Plot line widths from emission-line fits for DEEP2 sample
Created February 2019 but most recently added to MSC
"""

from chun_codes import systime, intersect

from astropy.io import ascii as asc
from astropy.io import fits
import astropy.units as u
import astropy.constants as const

from os.path import exists, join

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from glob import glob

from astropy.table import Table, vstack
from astropy import log

bbox_props = dict(boxstyle="square,pad=0.15", fc="white",
                  alpha=0.5, ec="none")


def exclude_outliers(objno):
    """
    Purpose
    --------
    This function excludes five specific sources that pass other exclusions
    Now included in MSC

    Parameters
    ----------
    objno -> the object id for the sources

    Returns
    -------
    flag -> list of ones and zeros that flags the bad sources
    """

    flag = np.zeros(len(objno), dtype=int)
    bad_data = np.array(['32007727', '32101412', '42006031', '32035286',
                         '14023705'])
    for ii in range(len(bad_data)):
        idx = [xx for xx in range(len(objno)) if bad_data[ii] == str(objno[xx])][0]
        print(ii, idx)
        flag[idx] = 1
    return flag


def get_data(out_path):
    """
    Purpose
    --------
    Used in the main function to collect the data needed 
    Parameters
    ----------
    out_path -> hardcoded into the main function (line 121) and is the path to get the datasets 

    Returns
    -------
    Data for main function 
    """

    path0 = '/Users/cly/data/DEEP2/DR4/f_current/'
    files = glob(path0+'*all_line_fit.fits')

    mass_file = '/Users/cly/Google Drive/MZEvolve/results_deeps.txt'
    mass_data = asc.read(mass_file)
    logM = np.log10(mass_data['best.stellar.m_star'].data)

    combine_stack_file = join(out_path, 'DEEP2_all_line_fit.fits')

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

    return data0, mass_data, logM
# enddef


def main(silent=False):

    """
    Main function to plot line widths of emission-line fits

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 25 February 2019
    """

    if not silent:
        log.info('### Begin main : '+systime())

    out_path = '/Users/cly/Google Drive/Zcalbase_gal/dataset/'

    data0, mass_data, logM = get_data(out_path)

    OIII = 1.33*data0['OIIIR_FLUX_MOD'].data
    OII = data0['OII_FLUX_MOD'].data
    HB = data0['HB_FLUX_MOD'].data

    SNR2_ini = data0['OII_SNR'].data
    SNR3_ini = data0['OIIIR_SNR'].data
    SNRH_ini = data0['HB_SNR'].data

    # exclude_flag = exclude_outliers(data0['OBJNO'])

    det3 = np.where((SNR2_ini >= 3) & (SNR3_ini >= 3) & (SNRH_ini >= 3) &
                    (OII > 0) & (OIII > 0) & (HB > 0))[0]
    print('# size det3 : ', len(det3))

    lR23 = np.log10((OII + OIII)/HB)
    lO32 = np.log10(OIII/OII)

    print('lR23 minmax : ', np.min(lR23[det3]), np.max(lR23[det3]))
    print('lO32 minmax : ', np.min(lO32[det3]), np.max(lO32[det3]))

    out_pdf = out_path + 'line_width.pdf'
    pp = PdfPages(out_pdf)

    # Velocity dispersion plot
    out_pdf1 = out_path + 'line_width_comparison.pdf'
    pp1 = PdfPages(out_pdf1)

    xlim = [0, 750.0]

    hist_bins = np.arange(xlim[0], xlim[1], 10)

    str_lines = ['OIIB', 'OIIR', 'HB', 'OIIIB', 'OIIIR']
    line_names = [r'[OII]$\lambda$3726', r'[OII]$\lambda$3728', r'H$\beta$',
                  r'[OIII]$\lambda$4959', r'[OIII]$\lambda$5007']
    lambda0 = [3726.16, 3728.91, 4861.32, 4958.91, 5006.84]

    c_value = const.c.to(u.km/u.s).value

    sig_limit = 250  # Dispersion limit in km/s
    # sig_arr   = np.zeros( (len(str_lines), len(data0)), dtype=np.int)
    sig_flag = np.zeros((len(str_lines), len(data0)), dtype=np.int)
    # sig_arr: contains flag for line detection
    # sig_flag: flags those with high dispersion

    # First loop to get those with high dispersion
    for ii in range(len(str_lines)):
        line = str_lines[ii]
        x_temp = data0[line+'_SIGMA'].data

        # Require S/N >= 3.0
        SN_cut = np.where(data0[line+'_SNR'].data >= 3.0)[0]

        x_temp = x_temp/lambda0[ii] * c_value

        t_cmd = line+'_vdisp = x_temp'
        exec(t_cmd)

        sig_idx = np.where(x_temp >= sig_limit)[0]
        sig_idx_SN = intersect(SN_cut, sig_idx)
        sig_flag[ii, sig_idx_SN] = 1

    sig_flag_total = np.sum(sig_flag, axis=0)

    high_disp = np.where(sig_flag_total >= 2)[0]
    print("high_disp : ", len(high_disp))

    # Plot velocity dispersion against each line
    fig_disp, ax_disp = plt.subplots(nrows=5, ncols=5)

    for ii in range(len(str_lines)):
        line = str_lines[ii]
        fig, ax = plt.subplots(nrows=2, ncols=2)

        x_temp = data0[line+'_SIGMA'].data

        good = np.where((x_temp < 90) & (x_temp > -90))[0]

        x_temp = x_temp/lambda0[ii] * c_value

        det3_good = intersect(det3, good)
        ax[0, 0].hist(x_temp[good], bins=hist_bins, alpha=0.5)

        ax[0, 0].hist(x_temp[det3_good], bins=hist_bins, color='red')

        ax[0, 0].annotate(line_names[ii], [0.95, 0.95], xycoords='axes fraction', ha='right', va='top')
        ax[0, 0].set_ylabel('N', fontsize=16)
        ax[0, 0].set_xticklabels([])

        ax[0, 1].scatter(x_temp[good], logM[good], alpha=0.5, edgecolor='none')
        ax[0, 1].scatter(x_temp[det3], logM[det3], facecolor='none', edgecolor='green', linewidth=0.5)

        ax[0, 1].scatter(x_temp[high_disp], logM[high_disp], marker='x', color='red')

        ax[0, 1].set_ylim([7.5, 12.0])
        ax[0, 1].set_xticklabels([])
        ax[0, 1].set_ylabel(r'$\log(M_{\star}/M_{\odot})$', fontsize=16)

        ax[1, 0].scatter(x_temp[good], lR23[good], alpha=0.5, edgecolor='none')
        ax[1, 0].scatter(x_temp[det3], lR23[det3], facecolor='none', edgecolor='green', linewidth=0.5)
        ax[1, 0].scatter(x_temp[high_disp], lR23[high_disp], marker='x', color='red')

        ax[1, 0].set_ylim(0.0, 2.25)
        ax[1, 0].set_xlabel(r'$\sigma$ [km s$^{-1}$]', fontsize=16)
        ax[1, 0].set_ylabel(r'$\log(R_{23})$', fontsize=16)

        ax[1, 1].scatter(x_temp[good], lO32[good], alpha=0.5, edgecolor='none')
        ax[1, 1].scatter(x_temp[det3], lO32[det3], facecolor='none', edgecolor='green', linewidth=0.5)
        ax[1, 1].scatter(x_temp[high_disp], lO32[high_disp], marker='x', color='red')

        ax[1, 1].set_ylim(-1, 2.25)
        ax[1, 1].set_xlabel(r'$\sigma$ [km s$^{-1}$]', fontsize=16)
        ax[1, 1].set_ylabel(r'$\log(O_{32})$', fontsize=16)

        for t_ax in ax.flatten():
            t_ax.set_xlim(xlim)

        fig.set_size_inches(11, 8.5)
        plt.subplots_adjust(left=0.06, right=0.99, top=0.97, bottom=0.065,
                            wspace=0.17, hspace=0.01)
        fig.savefig(pp, format='pdf')

        exec('vdisp_ii = '+line+'_vdisp')
        for jj in range(len(str_lines)):
            line_jj = str_lines[jj]
            if ii != jj:
                exec('vdisp_jj = '+line_jj+'_vdisp')
                in_range = np.where((vdisp_jj >= 0) & (vdisp_jj <= 2000) &
                                    (vdisp_ii >= 0) & (vdisp_ii <= 2000))[0]
                ax_disp[ii][jj].scatter(vdisp_jj[in_range], vdisp_ii[in_range], marker=',',
                                        s=1, color='k', edgecolor='none', alpha=0.5)

                pt2 = np.array([1, 1e3])
                ax_disp[ii][jj].plot(pt2, pt2, 'k--')

                t_diff = np.log10(vdisp_jj[in_range]/vdisp_ii[in_range])
                med0 = np.median(t_diff)
                avg0 = np.mean(t_diff)
                sig0 = np.std(t_diff)
                ax_disp[ii][jj].plot(pt2, pt2*10**(-med0), 'b-')
                ax_disp[ii][jj].plot(pt2, pt2*10**(-avg0), 'r-')
                ax_disp[ii][jj].plot(pt2, pt2*10**(-avg0-sig0), 'r--')
                ax_disp[ii][jj].plot(pt2, pt2*10**(-avg0+sig0), 'r--')

                a_txt = r'$\mu$'+' : %0.3f\n ' % avg0
                a_txt += r'$\tilde{\Delta}$'+' : %0.3f\n' % med0
                a_txt += r'$\sigma$ : %0.3f' % sig0
                ax_disp[ii][jj].annotate(a_txt, [0.975, 0.025], xycoords='axes fraction',
                                         ha='right', va='bottom', fontsize=6,
                                         bbox=bbox_props)

                ax_disp[ii][jj].set_xlim([25, 900])
                ax_disp[ii][jj].set_ylim([25, 900])

                ax_disp[ii][jj].set_yscale('log')
                ax_disp[ii][jj].set_xscale('log')

                if ii != len(str_lines)-1:
                    ax_disp[ii][jj].set_xticklabels([])
                if jj != 0:
                    ax_disp[ii][jj].set_yticklabels([])
            else:
                ax_disp[ii][jj].axis('off')
                ax_disp[0][jj].set_title(line_names[jj], fontsize=14)
                ax_disp[-1][jj].set_xlabel(r'$\sigma$ [km s$^{-1}$]')
                ax_disp[jj][0].set_ylabel(r'$\sigma$ [km s$^{-1}$]')
                ax_disp[jj][0].annotate(line_names[jj], [-0.5, 0.5], xycoords='axes fraction',
                                        ha='right', va='center', rotation=90, fontsize=14)
    # endfor

    log.info('Writing : '+out_pdf)
    pp.close()

    log.info('Writing : '+out_pdf1)
    fig_disp.set_size_inches(8.5, 8.5)
    fig_disp.subplots_adjust(left=0.12, right=0.995, top=0.96, bottom=0.065, wspace=0.01,
                             hspace=0.01)
    fig_disp.savefig(pp1, format='pdf')
    pp1.close()

    if not silent:
        log.info('### End main : '+systime())
