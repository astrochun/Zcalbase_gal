import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy.optimize import curve_fit
from os.path import join

from Metallicity_Stack_Commons.analysis.fitting import gauss
import general

fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'
O2_Histograms = join(fitspath, 'O2_Histograms.pdf')

# Spectral R23 and O32: Averages that are calculated from the flux calculations: spectral averages
spectral = '/Users/reagenleimbach/Desktop/Zcalbase_gal/April_15/combined_flux_table.tbl'
data1 = asc.read(spectral)

# Voronoi Outputs: R_23 and O_32 values for all galaxies with a column specifying bins
outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar21_voronoi_2d_binning_output.txt'
voronoi = np.genfromtxt(outfilevoronoi)


fluxes = [1e-17, 1e-16, 3e-17, 5e-16]
cont0 = [1e-18, 3e-18, 5e-18, 2e-17]
sig0 = [1.0, 2.3, 3.4, 2.1]


def test_gauss():
    x0 = np.arange(4950, 5050, 1)

    spec2d = np.zeros((len(fluxes), len(x0)))

    total_flux = np.zeros(len(fluxes))
    for ff in range(len(fluxes)):
        input = [5007, sig0[ff], fluxes[ff], cont0[ff]]
        spec2d[ff] = gauss(x0, *input)
        
        # o1,o2 = curve_fit(spec2d[ff], x0, fluxes[ff], p0=input)
        
        total_flux[ff] = np.sum(spec2d[ff] - cont0[ff])
        # fit_arr[ff] = fit
        plt.plot(x0, spec2d[ff])
        # plt.plot(x0, fit, color ='g', linestyle='dashed')

    print '## Avg fluxes : ', np.average(total_flux)

    stack0 = np.average(spec2d, axis=0)
    # print stack0.shape

    plt.plot(x0, stack0, color='magenta', linestyle='dashed')
    input = [5007, 3.0, np.average(fluxes), np.average(cont0)]
    o1, o2 = curve_fit(gauss, x0, stack0, p0=input)
    stack_fluxes = np.sum(gauss(x0, *o1) - o1[-1])
    print 'Composite flux : ', stack_fluxes


def single_bin_stacking():
    scalefact = 1e-17
    for ii in range(1, 5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])
    # fits.writeto(fitspath+'f3_0716/DEEP2_Fields_combined', data0)

    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']
     
    SNR2 = data0['OII_SNR']
    SNR3 = data0['OIIIR_SNR']
    SNRH = data0['HB_SNR']
    # SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                    (O2 > 0) & (O3 > 0) & (Hb>0))[0]

    data3 = data0[det3]
    O2 = data3['OII_FLUX_MOD']/scalefact
    O3 = 1.33*data3['OIIIR_FLUX_MOD']/scalefact
    Hb = data3['HB_FLUX_MOD']/scalefact

    O3_det = np.where(O3 >= (np.max(O3)-10))[0]
    print O3_det
    print len(O3_det)
    print np.max(O3)


def histogram_composite_plotting():
    scalefact = 1e-17
    binNum = voronoi[:, 4]
    
    O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3 = general.get_det3()
    binNum0 = list(set(binNum))

    OII = data1['OII_3727_Flux_Observed'].data
    OIII4959 = data1['OIII_4958_Flux_Observed'].data
    OIII5007 = data1['OIII_5007_Flux_Observed'].data
    H_BETA = data1['HBETA_Flux_Observed'].data
    N_Gal = data1['N_Galaxies'].data

    O2_arr = O2[det3]
    O3_arr = O3[det3]
    Hb_arr = Hb[det3]

    pdf_pages = PdfPages(O2_Histograms)
    for ii in range(len(binNum0)):
        # print N_Gal[ii]
        random_bin = binNum0[ii]
        bin_idx = [xx for xx in range(len(voronoi)) if binNum[xx] == random_bin]
        
        avg_OII = OII[ii]
        avg_OIII = 1.33 * OIII5007[ii]
        avg_Hb = H_BETA[ii]
        
        less_O2 = np.where((O2_arr[bin_idx] >= -90) & (O2_arr[bin_idx] <= 90))[0]
        O2_array = O2_arr[bin_idx]/scalefact
        not_O2 = N_Gal[ii] - len(less_O2)
        
        less_O3 = np.where((O3_arr[bin_idx] >= -90) & (O3_arr[bin_idx] <= 90))[0]
        O3_array = O3_arr[bin_idx]/scalefact
        not_O3 = N_Gal[ii] - len(less_O3)
        
        less_Hb = np.where((Hb_arr[bin_idx] >= -90) & (Hb_arr[bin_idx] <= 90))[0]
        Hb_array = Hb_arr[bin_idx]/scalefact
        not_Hb = N_Gal[ii] - len(less_Hb)
        
        # O2
        fig, ax_arr = plt.subplots()
        ax_arr.hist(O2_array, bins=20)
        ax_arr.axvline(x=avg_OII, color='r', linestyle='--', label='Composite Average')
        ax_arr.axvline(x=O2_array.mean(), color='g', linestyle='--', label='Histogram Average')
        ax_arr.set_ylabel('Number of Galaxies')
        ax_arr.set_xlabel('Fluxes (Scaled to 10^-17)')
        ax_arr.set_title('O2_bin' + str(ii))
        txt0 = r'N_Galaxies=%.3f' % (N_Gal[ii]) + '\n'
        txt0 += 'Number of Galaxies not used:%.3f' % not_O2
        ax_arr.annotate(txt0, [0.95, 0.80], xycoords='axes fraction', va='top', ha='right', fontsize='6')
        ax_arr.legend(loc=1,  fontsize='x-small')
        
        fig.savefig(pdf_pages, format='pdf')
        
        # O3
        fig, ax_arr = plt.subplots()
        ax_arr.hist(O3_array, bins=20)
        
        ax_arr.axvline(x=avg_OIII, color='r', linestyle='--', label='Composite Average')
        ax_arr.axvline(x=O3_array.mean(), color='g', linestyle='--', label='Histogram Average')
        ax_arr.set_ylabel('Number of Galaxies')
        ax_arr.set_xlabel('Fluxes (Scaled to 10^-17)')
        ax_arr.set_title('O3_bin' + str(ii))
        txt0 = r'N_Galaxies=%.3f' % (N_Gal[ii]) + '\n'
        txt0 += 'Number of Galaxies not used:%.3f' % not_O3
        ax_arr.annotate(txt0, [0.95, 0.80], xycoords='axes fraction', va='top', ha='right', fontsize='6')
        ax_arr.legend(loc=1,  fontsize='x-small')
        fig.savefig(pdf_pages, format='pdf')

        # Hb
        fig, ax_arr = plt.subplots()
        ax_arr.hist(Hb_array, bins=20)
        
        ax_arr.axvline(x=avg_Hb, color='r', linestyle='--', label='Composite Average')
        ax_arr.axvline(x=Hb_array.mean(), color='g', linestyle='--', label='Histogram Average')
        ax_arr.set_ylabel('Number of Galaxies')
        ax_arr.set_xlabel('Fluxes (Scaled to 10^-17)')
        ax_arr.set_title('Hb_bin' + str(ii))
        txt0 = r'N_Galaxies=%.3f' % (N_Gal[ii]) + '\n'
        txt0 += 'Number of Galaxies not used:%.3f' % not_Hb
        ax_arr.annotate(txt0, [0.95, 0.80], xycoords='axes fraction', va='top', ha='right', fontsize='6')
        ax_arr.legend(loc=1, fontsize='x-small')
        fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()
