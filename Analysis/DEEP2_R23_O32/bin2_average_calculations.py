#Checking voronoi bin2 and manually calculating R23 and O32 values

import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from os.path import exists
import numpy.ma as ma
from matplotlib.gridspec import GridSpec
from pylab import subplots_adjust
from astropy.convolution import Box1DKernel, convolve
from astropy.table import Table
from astropy.table import vstack 

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/'

tab= '/Users/reagenleimbach/Desktop/Zcalbase_gal/asc_table_voronoi.tbl'
asc_tab = asc.read(tab)



outfilevoronoi = '/Users/reagenleimbach/Desktop/Zcalbase_gal/Mar21_voronoi_2d_binning_output.txt'
voronoi = np.loadtxt(outfilevoronoi)
#voronoi = voronoi.transpose()
xcoor = [3726.16, 3728.91, 3797.90, 3835.38, 3868.74, 3889.05, 3888.65, 3967.51, 3970.07, 4340.46, 4363.21, 4471.5, 4958.91, 5006.84, 4101.73, 4363.21, 4861.32,5006.84]


RestframeMaster = r'/Users/reagenleimbach/Desktop/Zcalbase_gal/Master_Grid.fits'

def average_vals_bin2(): 
    for ii in range(1,5):
        file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
        data  = Table(fits.getdata(file1))
        if ii == 1:
            data0 = data
        else:
            data0 = vstack([data0, data])

        #print 'data0 : ', len(data0)
    O2 = data0['OII_FLUX_MOD']
    O3 = 1.33*data0['OIIIR_FLUX_MOD']
    Hb = data0['HB_FLUX_MOD']
   
    R23 = (O2+O3)/Hb
    O32 = O3/O2
    #print R23, O32
    
    SNR2 = data0['OII_SNR']
    SNR3 = data0['OIIIR_SNR']
    SNRH = data0['HB_SNR']

    #SNR code: This rules out major outliers by only using specified data
    det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                    (O2 > 0) & (O3 > 0) & (Hb>0))[0]
    print len(det3)

    data3 = data0[det3]
    O2 = data3['OII_FLUX_MOD']
    O3 = 1.33*data3['OIIIR_FLUX_MOD']
    Hb = data3['HB_FLUX_MOD']

    lR23 = voronoi[:,2]
    lO32 = voronoi[:,3]
    binnum = voronoi[:,4]
    binnum0= list(set(binnum))
    
    avg_l_R23 = np.zeros(len(binnum0))
    avg_l_O32 = np.zeros(len(binnum0))
    for ii in binnum0:
        idx = [xx for xx in range(len(det3)) if binnum[xx] == ii]
            
        avg_O2 = np.average(O2[idx])
        avg_O3 = np.average(O3[idx])
        avg_Hb = np.average(Hb[idx])

        avg_l_R23[ii] = np.log10((avg_O2 +avg_O3)/avg_Hb)
        avg_l_O32[ii] = np.log10(avg_O3/avg_O2)

    return avg_l_R23, avg_l_O32



