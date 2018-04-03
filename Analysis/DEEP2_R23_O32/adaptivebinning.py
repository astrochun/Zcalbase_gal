##Adaptive binning for DEEP2 data
##Running Voronoi code 

'''
   Notes
   -----
   Modified by Chun Ly, 15 September 2017
    - Call voronoi_2d_binning with currentBin input
   Modified by Reagen Leimbach, 5 October 2017
    - Changed location of vornoi_2d_binning with relationship to adaptivebinning 
      and commented out all unnecessary print statments to fix bug
    - Add Ascii table that has the values calculated by voronoi_2d_binning 

'''

from os import path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
import scipy 
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.table import vstack
from pprint import pprint


import voronoi_2d_binning

fitspath='/Users/reagenleimbach/Desktop/Zcalbase_gal/' 

#Creates a fits table 
for ii in range(1,5):
    file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
    data  = Table(fits.getdata(file1))
    if ii == 1:
        data0 = data
    else:
        data0 = vstack([data0, data])


O2 = data0['OII_FLUX_MOD'].data
O3 = 1.33*data0['OIIIR_FLUX_MOD'].data
Hb = data0['HB_FLUX_MOD'].data
R23 = (O2+O3)/Hb
O32 = (O3/O2)

O2_Noise = data0['OII_NOISE']

SNR2 = data0['OII_SNR']
SNR3 = data0['OIIIR_SNR']
SNRH = data0['HB_SNR']

#SNR code: This rules out major outliers by only using specified data
det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                (O2 > 0) & (O3 > 0) & (Hb>0))[0]

#print np.min(R23[det3]), np.max(R23[det3]), np.mean(R23[det3]), np.min(O32[det3]), np.max(O32[det3]), np.mean(O32[det3])
#print len(det3)



def voronoi_binning_DEEP2():
    #file_dir = path.dirname(path.realpath(__file__))  # path of this procedure, edit this to correct path

    #x, y, signal, noise = (fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits', upack=1, skiprows=3)
    #targetSN= 1.00
    #what would these signal and noise values be? 

    # Perform the actual computation. The vectors
    # (binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale) #Can output any of these 
    # are all generated in *output*

    print voronoi_2d_binning.__file__
    
    #print len(R23),len(O32)
    signal = np.repeat(1, len(det3))
    noise = np.repeat(1, len(det3))
    #print np.min(R23[det3]), np.max(R23[det3])

    # + on 15/09/2017
    rR23 = R23[det3]
    rO32 = O32[det3]
    lR23 = np.log10(R23[det3])
    lO32 = np.log10(O32[det3])
    avg_R23 = np.average(lR23)
    avg_O32 = np.average(lO32)

    dist0      = np.sqrt((lR23-avg_R23)**2 + (lO32-avg_O32)**2)
    currentBin = np.argmin(dist0)
    print currentBin

    #binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale
    test = voronoi_2d_binning.voronoi_2d_binning(lR23, lO32, signal, noise, 10, plot=True, pixelsize=0.1,quiet=False, currentBin=currentBin)
    #print(test.shape)

    print np.sum(test[6])
    # Save to a text file the initial coordinates of each pixel together
    # with the corresponding bin number computed by this procedure.
    # binNum uniquely specifies the bins and for this reason it is the only
    # number required for any subsequent calculation on the bins.
    #
    np.savetxt(fitspath+'Mar21_voronoi_2d_binning_output.txt', np.column_stack([rR23,rO32, lR23, lO32, test[0]])) #test[0]=classe = binnumber 
    #fmt='%6.3f %6.3f %5i')
    
    #n=('xnode','ynode','xBar','yBar','sn','area','scale')
    #tab0= Table([test[1], test[2], test[3], test[4], test[5], test[6], test[7]], names=n)
    #asc.write(tab0, fitspath+'asc_table_voronoi.tbl')
    #print tab0

    #return test
#-----------------------------------------------------------------------------

#if __name__ == '__main__':
#
#    voronoi_binning_DEEP2()



