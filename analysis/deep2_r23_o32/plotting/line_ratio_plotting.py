import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join


def plotting_data1(fitspath, dataset, combine_flux_file, bin_info_file):
    """
    Plots measurements to check to make sure that data follows
    a one-to-one line

    :param fitspath: str. path where files are called from and saved to
    :param dataset: str.  keyword used to define binning method
    :param combine_flux_file: str. table of data calculated for emission line
    :param bin_info_file: str. table of averages calculated in binning code

    PDF File: fitspath + 'line_ratio_plots.pdf'
    No returns
    """

    pdf_file = join(fitspath, 'line_ratio_plots.pdf')
    pp = PdfPages(pdf_file)

    print("### combine_flux_file : " + combine_flux_file)
    fitted_data = asc.read(combine_flux_file)

    print("### bin_info_file : " + bin_info_file)
    bin_info_tab = asc.read(bin_info_file)
    OII = fitted_data['OII_3727_Flux_Observed']
    OIII5007 = fitted_data['OIII_5007_Flux_Observed']
    H_BETA = fitted_data['HBETA_Flux_Observed']
    binnum = fitted_data['N_stack']
    ID = fitted_data['bin_ID']
    print('binnum:', binnum, len(binnum))

    R23_composite = np.zeros(binnum.shape[0])
    O32_composite = np.zeros(binnum.shape[0]) 
    for ii in range(len(binnum)):
        R23_comp = np.log10((OII[ii] + (1.33 * OIII5007[ii]))/H_BETA[ii])
        O32_comp = np.log10((1.33 * OIII5007[ii])/OII[ii])
        print(R23_comp, O32_comp)
        R23_composite[ii] = R23_comp
        O32_composite[ii] = O32_comp
    
    R23_raw = bin_info_tab['logR23_avg']
    O32_raw = bin_info_tab['logO32_avg']
    binnum_raw = bin_info_tab['N_stack']

    if dataset != 'Grid':
        for rr in range(len(binnum)):
            if binnum[rr] == binnum_raw[rr]:
                print('equal', binnum[rr], binnum_raw[rr])

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(R23_raw, R23_composite, marker='o', facecolor='none',
                   edgecolor='b', label='R23 Ratio: Voronoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset + ' Raw vs. Composite for R23')
    for rr in range(len(ID)):
        ax_arr.annotate(ID[rr], (R23_raw[rr], R23_composite[rr]))
    ax_arr.set_xlabel(r'Raw log($R_{23}$)')
    ax_arr.set_ylabel(r'Composite log($R_{23}$)')
    ax_arr.plot([0.0, 1.3], [0.0, 1.3], 'k-')
    
    fig.savefig(pp, format='pdf')

    fig, ax_arr = plt.subplots()
    ax_arr.scatter(O32_raw, O32_composite, marker='o', facecolor='none',
                   edgecolor='b', label='O32 Ratio: Voronoi Raw vs. Composite')
    ax_arr.legend(loc=0)
    ax_arr.set_title(dataset + 'Raw vs. Composite for O32')
    for oo in range(len(ID)):
        ax_arr.annotate(ID[oo], (O32_raw[oo], O32_composite[oo]))
    ax_arr.set_xlabel(r'Raw log($O_{32}$)')
    ax_arr.set_ylabel(r'Composite log($O_{32}$)')

    ax_arr.plot([-1, 1.2], [-1, 1.2], 'k-')
    fig.savefig(pp, format='pdf')

    pp.close()

    fig.clear()
