
# NOT CALLED IN GENERAL FUNCTIONS
# MAKES PLOTS


# Creating a color bar for the R23 and O32 values

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages

# Voronoi 14
'''fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi14_1016/'
asc_table1 = fitspath+'voronoi14_binning_averages.tbl'
asc_table2 = fitspath+ 'voronoi14_2d_binning_datadet3.tbl'

temp_tab = fitspath + 'Voronoi14_temperatures_metalicity.tbl'''

# Voronoi 10
'''fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/Voronoi10_1018/'
asc_table1 = fitspath+'voronoi10_binning_averages.tbl'
asc_table2 = fitspath+ 'voronoi10_2d_binning_datadet3.tbl'

temp_tab = fitspath + 'Voronoi10_temperatures_metalicity.tbl'''

# Grid
fitspath = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Gridmethod_1018/'
# asc_table1 = fitspath + 'voronoi14_binning_averages.tbl'
asc_table1 = fitspath + 'Gridbinning_averages.tbl'

temp_tab = fitspath + 'Grid_temperatures_metalicity.tbl'

asc_table1 = asc.read(asc_table1)
# asc_table = asc.read(asc_table2)
temp_table = asc.read(temp_tab)


def color_for_bin():
    targetSN = 14
    xBar = asc_table['xBar']
    yBar = asc_table['yBar']
    xnode = asc_table['xnode']
    ynode = asc_table['ynode']
    area = asc_table['area']

    w = area == 1
    plt.clf()
    plt.subplot(211)
    rnd = np.argsort(np.random.random(xnode.size))  # Randomize bin colors
    #_display_pixels(x, y, rnd[classe], pixelsize)
    # plt.plot(xnode, ynode, '+w', scalex=False, scaley=False) # do not rescale after imshow()
    # plt.xlabel('R (arcsec)')
    # plt.ylabel('R (arcsec)')
    # plt.title('Map of Voronoi bins')
    
    plt.subplot(212)
    rad = np.sqrt(xBar**2 + yBar**2)  # Use centroids, NOT generators
    plt.plot(rad[~w], sn[~w], 'or', label='Voronoi bins')
    plt.xlabel('R (arcsec)')
    plt.ylabel('Bin S/N')
    plt.axis([np.min(rad), np.max(rad), 0, np.max(sn)])  # x0, x1, y0, y1
    if np.sum(w) > 0:
        plt.plot(rad[w], sn[w], 'xb', label='single spaxels')
    plt.axhline(targetSN)
    plt.legend()
    plt.pause(1) 


def R23vsO32_plot():
    pdf_name = 'R23vsO32_color_comandavg.pdf'
    pdf_pages = PdfPages(fitspath+pdf_name)

    # logR23 = asc_table['log(R23)']
    # logO32 = asc_table['log(O32)']
    com_O_log= temp_table['com_O_log']
    temp = temp_table['Temperature']
    R23_com = temp_table['R23_Composite']
    O32_com = temp_table['O32_Composite']
    n_gal = temp_table['N_Galaxies']
    
    xBar = asc_table1['xBar']
    yBar = asc_table1['yBar']

    fig1, ax1 = plt.subplots()
    cm = plt.cm.get_cmap('BuPu_r')

    vmin = np.nanmin(com_O_log)  # temp)
    vmax = np.nanmax(com_O_log[com_O_log != np.inf])   # temp
    c = com_O_log  # temp

    print('min, max, c: ', vmin, vmax, c)
    # scat= plt.scatter(logR23, logO32, c=logR23, vmin=min(logR23), vmax=max(logR23), marker='.', cmap=cm)
    scat = plt.scatter(xBar, yBar, c=c, vmin=vmin, vmax=vmax, marker='.', cmap=cm)
    scat_com = plt.scatter(R23_com, O32_com, c=c, vmin=vmin, vmax=vmax, marker='*', cmap=cm)
    cax = plt.axes([0.9, 0.15, 0.8, 0.75])
    cbar = plt.colorbar(scat, cax=cax, orientation='vertical', label='Oxygen Abundance')
    for rr in range(len(n_gal)):
        ax1.annotate(str(n_gal[rr]), (xBar[rr], yBar[rr]), xycoords='data', fontsize=5)
    for oo in range(len(n_gal)):
        ax1.annotate(str(n_gal[oo]), (R23_com[oo], O32_com[oo]), xycoords='data', fontsize=5, fontweight='bold')
    ax1.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax1.set_xlabel(r'log($R_{23}$)')
    ax1.set_ylabel(r'log($O_{32}$)')

    #fig1.set_size_inches(8,8)

    pdf_pages.savefig()
    pdf_pages.close()
