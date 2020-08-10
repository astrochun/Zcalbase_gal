import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join

from Metallicity_Stack_Commons.column_names import filename_dict


def graph_bins(fitspath, n_split, bin_info, individual_info, pdf_file):
    """
    Purpose: This function takes inputed grid parameters and graphs them on the log(R23)--log(O32) plane 
    Inputs:
           n_split    -> from general.py input used for calculating bins
           ascii_file -> file with the minimum, maximum, and average R23 and O32 value
           pdf_file   -> name of pdf file made by code

    Returns:
           none
    Outputs:
        pdf_file with the saved figure (possible naming convention 'manual_binning_shading.pdf'

    Currently works for manual binning with an n_split (see general function) of 3. More updates would
    need to be made for a different n_split to correct indexing through arrays. 
    """

    bin_info = join(fitspath, filename_dict['bin_info'])
    individual_info = join(fitspath, filename_dict['individual_prop'])
    pdf_file = join(fitspath, 'manual_binning_shading.pdf')

    pdf_pages = PdfPages(pdf_file)

    tab = asc.read(bin_info)
    idv_tab = asc.read(individual_info)

    lR23_min = tab['logR23_min'].data
    lO32_min = tab['logO32_min'].data
    lR23_avg = tab['logR23_avg'].data
    lO32_avg = tab['logO32_avg'].data
    lR23_max = tab['logR23_max'].data
    lO32_max = tab['logO32_max'].data

    lR23 = idv_tab['logR23'].data
    lO32 = idv_tab['logO32'].data
    
    fig, ax = plt.subplots()

    count = 0
    for aa in range(len(lR23_min)):
        xmin = lR23_min[aa]
        if aa <= (len(lO32_min) - n_split - 1):
            xmax = lR23_max[aa]
        else:
            xmax = lR23_max[-1]
        
        x_value = [xmin, xmax]
        y_average = [lO32_avg[aa], lO32_avg[aa]]

        y_value = [lO32_min[aa], lO32_max[aa]]
        x_average = [lR23_avg[aa], lR23_avg[aa]]
        
        plt.plot(x_average, y_value, '--', linewidth=0.75, color='k')
        plt.plot(x_value, y_average, '--', linewidth=0.75, color='k')

        # for cc in range(3):
        if aa <= (len(lR23_min) - 2):
            xx = [xmin, lR23_max[aa]]
        else:
            xx = [lR23_min[aa-1], lR23_max[-1]]
        if aa <= (len(lO32_max) - 1):

            ymax = lO32_max[aa]
        else:
            ymax = lO32_max[-1]
        print('ymin, ymax = ', lO32_min[aa], ymax)
        ax.fill_between(xx, lO32_min[aa], ymax, alpha=0.35)
        count += 1

    print(count)
    
    finite0 = np.where((np.isfinite(lR23)) & (np.isfinite(lO32)))[0]
    x1 = lR23[finite0]
    y1 = lO32[finite0]
    ax.scatter(x1, y1, 0.75, facecolor='b', edgecolor='face', marker='*', alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for Manual Binning')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()
