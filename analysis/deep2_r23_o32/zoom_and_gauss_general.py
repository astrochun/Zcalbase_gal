"""
PERFORMS THE GAUSSIAN FITTING ON ALL EMISSIONS LINES
CALLED EVERY TIME THE GENERAL FUNCTIONS ARE RUN

Debugging Note:
   If 'x0' is infeasible error occurs, check the para_bound values to
   make sure the expected values are within the range set up upper and
   lower limits.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column, hstack
from pylab import subplots_adjust
from scipy.optimize import curve_fit
from os.path import join

from Metallicity_Stack_Commons.analysis.fitting import gauss, double_gauss, oxy2_gauss
from Metallicity_Stack_Commons.analysis.fitting import movingaverage_box1D, rms_func
from Metallicity_Stack_Commons import lambda0, line_name, line_type
from Metallicity_Stack_Commons.column_names import filename_dict

from . import name_dict
from .log_commons import log_stdout

con1 = 3728.91 / 3726.16


def line_flag_check(fitspath, working_wave, lineflag, wave, y_norm,
                    line_name0, row, col, fig, ax_arr, log=None):
    """
    Purpose
      Plots a zoomed in plot of emission lines to check visually if the
      emission lines are excluded in the line flag.
      Note: Call is currently commented out
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    # New plots for line flagging
    out_pdf = join(fitspath, f"_lineflag_check_{line_name0}.pdf")
    pdf_pages2 = PdfPages(out_pdf)

    t_ax2 = ax_arr[row, col]
    t_ax2.plot(wave, y_norm, 'k', linewidth=0.4, label='Emission')
    t_ax2.set_xlim([working_wave + 150, working_wave - 45])
    t_ax2.plot(wave, lineflag, 'g', linestyle='dashed', label='Lineflag')

    fig.set_size_inches(8, 8)
    log.debug(f"Writing: {out_pdf}")
    fig.savefig(pdf_pages2, format='pdf')
    pdf_pages2.close()

    log.debug("finished.")


def get_gaussian_fit(dataset, s2, working_wave, x0, y_norm, x_idx, rms,
                     line_type0, log=None):
    """
    Purpose:
      Calculates the gaussian to fit each emission line using curve_fit

    Debugging Note:
    If 'x0' is infeasible error occurs, check the para_bound values to
    make sure the expected values are within the range set up upper and
    lower limits.
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    med0 = np.median(y_norm[x_idx])
    max0 = np.max(y_norm[x_idx]) - med0
    sigma = np.repeat(rms, len(x0))

    fail = 0
    o1 = None

    # curve_fit(gauss function, x0 ,y0, initial_guesses, sigma, parameter_bounds)
    # Single Emission Line
    if line_type0 == 'Single':
        p0 = [working_wave, 1.0, max0, med0]  # must have some reasonable values

        para_bound = ((working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0)),
                      (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0)))

        try:
            o1, o2 = curve_fit(gauss, x0[x_idx], y_norm[x_idx], p0=p0,
                               sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            log.warning('fail')
            fail = 1

    # Double Balmer Emission Line
    if line_type0 == 'Balmer':
        para_bound = (working_wave - 3.0, 0.0, 0.0, med0 - 0.05 * np.abs(med0), 0.0, -med0), \
                     (working_wave + 3.0, 10.0, 100.0, med0 + 0.05 * np.abs(med0), 30.0, 0.0)
        if dataset in ['R23_Grid', 'O32_Grid']:
            # must have some reasonable values
            p0 = [working_wave, 1.0, max0, med0, s2, -0.05*max0]

        if dataset == 'Grid':
            # must have some reasonable values
            p0 = [working_wave, 1.0, max0, med0, s2, -0.25*med0]

        if dataset in ['Voronoi10', 'Voronoi14', 'Voronoi20',
                       'Double_Bin', 'n_Bins']:
            # must have some reasonable values
            p0 = [working_wave, 1.0, max0, med0, s2, -0.5*med0]

        # Attempt fit
        try:
            o1, o2 = curve_fit(double_gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            log.warning('fail')
            fail = 1

    # OxygenII Emission Line
    if line_type0 == 'Oxy2':
        p0 = [working_wave, 1.0, 0.75*max0, med0, 1.0, max0]  # must have some reasonable values
        para_bound = (working_wave-3.0, 0.0, 0.0, med0-0.05*np.abs(med0), 0.0, 0.0),\
                     (working_wave+3.0, 10.0, 100.0, med0+0.05*np.abs(med0), 10.0, 100.0)

        try:
            o1, o2 = curve_fit(oxy2_gauss, x0[x_idx], y_norm[x_idx],
                               p0=p0, sigma=sigma[x_idx], bounds=para_bound)
        except ValueError:
            log.warning('fail')
            fail = 1

    log.debug("finished.")

    return o1, med0, max0


def equi_width_func(pos_comp, neg0, gauss0, x0, wave, y_norm):
    """
    Purpose:
      Equivalent width correction/computation


    Notes:
      bottom side of the iceburg / continuum
      take negative component of gauss and subtract off the positive component

      total gauss - the positive component?
      x = double gaussian of o1[0,1,2] = 0
      x-o1[3] from [-2.5sigma to 2.5 sigma]
      equivalent width in terms of Angstroms
      update plots
    """

    plt.plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
    plt.plot(x0, neg0, 'b', linewidth=0.25, label='Negative Component')
    plt.plot(x0, pos_comp, 'r', linewidth=0.25, label='Positive Component')
    plt.plot(x0, gauss0, 'g', linewidth=0.25, label='Gauss Fit')


# For each individual stack
# Electron temperature and the R23 and O32 values
# Plotting Zoomed
def zoom_gauss_plot(dataset, stack2d, dispersion, s2, wave,
                    working_wave, lineflag, bin_info_tab,
                    y_correction='', line_type='',
                    out_pdf='', line_name='', log=None):
    """
    Purpose
    ----------
    Main function that is called by run function (zm_general). Gets the data, fits a gaussian curve
    to each emission line, plots emission lines for each bin, and saves off .tbl files with curve information.

    Debugging Note:
    If 'x0' is infeasible error occurs, check the para_bound values to
    make sure the expected values are within the range set up upper and
    lower limits.
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    pdf_pages = PdfPages(out_pdf)

    nrows = 4
    ncols = 4
    x_idx = np.where((wave >= (working_wave-100)) &
                     (wave <= (working_wave+100)))[0]
    x0 = wave
    scalefact = 1e-17

    bin_id = bin_info_tab['bin_ID'].data

    # Initializing Arrays
    flux_g_array = np.zeros(stack2d.shape[0])
    flux_s_array = np.zeros(stack2d.shape[0])
    flux_neg_array = np.zeros(stack2d.shape[0])
    sigma1_array = np.zeros(stack2d.shape[0])
    median_array = np.zeros(stack2d.shape[0])
    norm_array = np.zeros(stack2d.shape[0])
    rms_array = np.zeros(stack2d.shape[0])
    SN_array = np.zeros(stack2d.shape[0])

    # Initializing Arrays for Balmer Graphing
    xbar_array = np.zeros(stack2d.shape[0])
    pos_amp_array = np.zeros(stack2d.shape[0])
    sig2_array = np.zeros(stack2d.shape[0])
    neg_amp_array = np.zeros(stack2d.shape[0])

    for rr in range(stack2d.shape[0]):
        y0 = stack2d[rr]
        y_norm = y0/scalefact

        row = rr // nrows % ncols
        col = rr % ncols
        if rr % (nrows*ncols) == 0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols,
                                       squeeze=False)

        t_ax = ax_arr[row, col]

        x1 = working_wave-100
        x2 = working_wave+100

        y_smooth = movingaverage_box1D(y_norm, 4, boundary='extend')

        rms_ang = rms_func(wave, dispersion, working_wave, y0, 0, lineflag)

        if y_correction:
            o1, med0, max0 = get_gaussian_fit(dataset, s2, working_wave, x0,
                                              y_smooth, x_idx, rms_ang, line_type,
                                              log=log)
        else:
            o1, med0, max0 = get_gaussian_fit(dataset, s2, working_wave, x0,
                                              y_norm, x_idx, rms_ang, line_type,
                                              log=log)

        # Calculating Flux: Signal Line Fit
        if o1 is not None:
            dx = x0[2]-x0[1]
            if line_type == 'Single':
                x_sigsnip = np.where((np.abs((x0 - working_wave))/o1[1]) <= 2.5)[0]
                gauss0 = gauss(x0, *o1)

            if line_type == 'Balmer':
                x_sigsnip = np.where(np.abs((x0 - working_wave))/o1[1] <= 2.5)[0]
                gauss0 = double_gauss(x0, *o1)

                o1_neg = [o1[0], o1[4], o1[5], o1[3]]
                neg0 = gauss(x0, *o1_neg)
                gauss0_diff = gauss0 - neg0
                y_norm_diff = y_norm[x_sigsnip] - neg0[x_sigsnip]

                # Equivalent Width calculation
                pos_comp = gauss0 - neg0
                neg_comp = neg0[x_sigsnip]
                flux_neg = np.sum(neg_comp*dx)
                equiv_wid = flux_neg/o1[3]

                # Keeping this here for possible later use. Might belong in MSC balmer module
                # pdfpages3 = PdfPages(fitspath + '/' + dataset + '_PlottingBalmer_' + line_name+'.pdf')
                # equi_width_func(fitspath, dataset, working_wave, pos_comp, neg0, gauss0, x0, wave,
                # y_norm, line_type, pdfpages3)

            if line_type == 'Oxy2':
                x_sigsnip = np.where(((x0-working_wave)/o1[1] >= -2.5) &
                                     ((x0-working_wave*con1)/o1[4] <= 2.5))[0]
                gauss0 = oxy2_gauss(x0, *o1)

            # Get fluxes
            if line_type in ['Single', 'Oxy2']:
                # Flux from gaussian distribution
                flux_g = np.sum((gauss0-o1[3])*dx)
                # Flux from snipping method (spectral flux)where snip off sigma > 2.5
                flux_s = np.sum((y_norm[x_sigsnip] - o1[3]) * dx)

            if line_type == 'Balmer':
                flux_g = np.sum(gauss0_diff * dx)
                flux_s = np.sum(y_norm_diff * dx)

            # Calculating rms
            rms_tot, rms_pix = rms_func(wave, dispersion, working_wave, y0, o1[1],
                                        lineflag)

            # Line Flag Checking Plots
            # line_flag_check(fitspath, working_wave, lineflag, wave, y_norm, stack2d,
            #                line_name, row, col, fig, ax_arr, log=log)

            # Array Population
            norm_array[rr] = max0
            rms_array[rr] = rms_tot
            SN_array[rr] = (flux_s/rms_tot)
            if line_type == 'Balmer':
                flux_neg_array[rr] = flux_neg

            flux_g_array[rr] = flux_g
            flux_s_array[rr] = flux_s
            xbar_array[rr] = o1[0]   # referred to as Center as well
            sigma1_array[rr] = o1[1]
            pos_amp_array[rr] = o1[2]
            median_array[rr] = o1[3]
            norm_array[rr] = max0
            rms_array[rr] = rms_tot
            SN_array[rr] = (flux_s/rms_tot)

            if line_type in ['Balmer', 'Oxy2']:
                sig2_array[rr] = o1[4]
                neg_amp_array[rr] = o1[5]

            # Residuals
            if line_type == 'Oxy2':
                x_sigsnip_2 = np.where(np.abs(x0 - 3727)/o1[1] <= 4.0)[0]
                resid = y_norm[x_sigsnip_2] - gauss0[x_sigsnip_2] + o1[3]
            else:
                x_sigsnip_2 = np.where(np.abs((x0 - working_wave))/o1[1] <= 3.0)[0]
                resid = y_norm[x_sigsnip_2] - gauss0[x_sigsnip_2] + o1[3]

            # Plotting
            if y_correction:
                t_ax.plot(wave, y_smooth, 'k', linewidth=0.3, label='Emission')
                t_ax.plot(x0, gauss0, 'm', linewidth=0.25, label='Gauss Fit')
            else:
                t_ax.plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
                t_ax.plot(x0, gauss0, 'b', linewidth=0.25, label='Gauss Fit')

            t_ax.plot(x0[x_sigsnip_2], resid, 'r', linestyle='dashed',
                      linewidth=0.2, label='Residuals')
            t_ax.set_xlim(x1 + 50, x2 - 50)

            if line_name == 'OIII_4363':
                t_ax.set_ylim(0, 1)

            if dataset in ['Grid', 'O32_Grid', 'R23_Grid', 'Double_Bin', 'n_Bins']:
                txt0 = 'Line: %.3f, ID: %i, R_23: %.3f O_32: %.3f\n' % \
                       (o1[0], bin_id[rr], bin_info_tab['logR23_min'][rr],
                        bin_info_tab['logO32_min'][rr]) + '\n'
                txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %i\n' % (rms_tot, rms_pix,
                                                              bin_info_tab['N_stack'][rr])
                if line_type == 'Balmer':
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f' % (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f' % (o1[2], o1[4], o1[5]) + '\n'
                if line_type in ['Single', 'Oxy2']:
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f' % (o1[3], o1[1], max0, o1[2]) + '\n'

                txt0 += 'F_G: %.3f F_S: %.3f' % (flux_g, flux_s) + '\n'
                txt0 += 'S/N: %.3f' % (SN_array[rr])
            else:
                txt0 = r'Line: %.3f, ID: %i  xnode=%.3f  ynode=%.3f' % \
                       (o1[0], bin_id[rr], bin_info_tab['logR23_min'][rr],
                        bin_info_tab['logO32_min'][rr]) + '\n'
                txt0 += 'R_23: %.3f O_32: %.3f\n' % (bin_info_tab['logR23_avg'][rr],
                                                     bin_info_tab['logO32_avg'][rr])
                txt0 += 'RMS: %.3f RMS/pix: %.3f, N: %i \n' % (rms_tot, rms_pix, bin_info_tab['N_stack'][rr])
                if line_type == 'Balmer':
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f' % (o1[3], o1[1], max0) + '\n'
                    txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f' % (o1[2], o1[4], o1[5]) + '\n'
                if line_type in ['Single', 'Oxy2']:
                    txt0 += r'Median: %.3f $\sigma$: %.3f  Norm: %.3f o1[2]: %.3f' % \
                            (o1[3], o1[1], max0, o1[2]) + '\n'
                txt0 += 'Flux_G: %.3f Flux_S: %.3f' % (flux_g, flux_s) + '\n'
                txt0 += 'S/N: %.3f' % (SN_array[rr])

            t_ax.annotate(txt0, [0.95, 0.95], xycoords='axes fraction',
                          va='top', ha='right', fontsize='5')

            for x in lambda0:
                t_ax.axvline(x=x, linewidth=0.15, color='k', linestyle='--')

            txt1 = 'Intensity'  # Units: r'$ergs *s^{-1} *cm^{-2} *\AA$'

            if col == 0:
                t_ax.set_ylabel(txt1)
            else:
                t_ax.set_yticklabels([])  # sets y-tick labels

            if row != nrows-1 and rr != stack2d.shape[0]-1:
                t_ax.set_xticklabels([])
        else:
            t_ax.plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
            t_ax.set_xlim(x1 + 50, x2 - 50)

        if (rr % (nrows * ncols) == nrows * ncols - 1) or rr == stack2d.shape[0] - 1:
            subplots_adjust(left=0.1, right=0.98, bottom=0.06, top=0.97,
                            hspace=0.05)

            fig.set_size_inches(8, 8)
            fig.savefig(pdf_pages, format='pdf')
            fig.clear()

    # Repetitive Columns have been added: The last four to six columns will be used for individual graphing

    if line_type == 'Single':
        n = ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS',
             'S/N', 'Center', 'Pos_Amp')

        n = tuple([f"{line_name}_{val}" for val in n])
        tab0 = Table([flux_g_array, flux_s_array, sigma1_array, median_array,
                      norm_array, rms_array, SN_array, xbar_array, pos_amp_array],
                     names=n)

    if line_type in ['Balmer', 'Oxy2']:
        n = ('Flux_Gaussian', 'Flux_Observed', 'Sigma', 'Median', 'Norm', 'RMS',
              'S/N', 'Center', 'Pos_Amp', 'Abs_Sigma', 'Abs_Norm')
        n = tuple([f"{line_name}_{val}" for val in n])

        tab0 = Table([flux_g_array, flux_s_array, sigma1_array, median_array,
                      norm_array, rms_array, SN_array, xbar_array,
                      pos_amp_array,  sig2_array, neg_amp_array], names=n)

        if line_type == 'Balmer':
            log.info('Adding an Equ_Width Column')
            names = f"EW_{np.int(working_wave)}_abs"
            equ_add = Column(name=names, data=flux_neg_array)
            tab0.add_column(equ_add, 2)

    pdf_pages.close()
    # pdfpages3.close() # Add back in equivalent width plots are added
    fig.clear()

    log.debug("finished.")
    return tab0


def zm_general(dataset, fitspath, stack2d, wave, lineflag, dispersion,
               y_correction, log=None):
    """
    Purpose
    ----------
    Run function for gaussian fitting step
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    s2 = 5.0  # a fitting requirement

    # Importing Bins Values
    bin_info_file = join(fitspath, filename_dict['bin_info'])
    log.info(f"Reading: {bin_info_file}")
    bin_info_tab = asc.read(bin_info_file)

    n2 = ('bin_ID', 'logR23_avg', 'logO32_avg', 'N_stack')
    table_stack = bin_info_tab[n2]

    for ii in range(len(lambda0)):
        out_pdf = join(fitspath, f"Zoomed_Gauss_{line_name[ii]}.pdf")
        log.info(f"out_pdf: {out_pdf}")

        em_tab = zoom_gauss_plot(dataset, stack2d, dispersion, s2,
                                 wave, lambda0[ii], lineflag, bin_info_tab,
                                 y_correction=y_correction, line_type=line_type[ii],
                                 out_pdf=out_pdf, line_name=line_name[ii],
                                 log=log)

        table_stack = hstack([table_stack, em_tab])

    out_ascii = join(fitspath, filename_dict['bin_fit'])
    log.info(f"Writing: {out_ascii}")
    asc.write(table_stack, out_ascii, format='fixed_width_two_line',
              overwrite=True)

    log.debug("finished.")
