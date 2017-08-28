#Currently binning with bin=0.25
#Numpy binary file or FITS binary table: read up on both and attempt both
#Change name of pdf, name of output file, bin size

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.table import vstack 

fitspath='/astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/' 
#pdf_pages = PdfPages(fitspath+'R23_O32_bin01_scatter_and_hexbin_MasterGrid.pdf') #open pdf document
pdf_pages = PdfPages('astrochun/Zcalbase_gal/Analysis/DEEP2_R23_O32/R23_O32_bin025_scatter_and_hexbin.pdf'
#marker0 = ['b','g','r','m']

#Creates a fits table 
for ii in range(1,5):
    file1 = fitspath+'f3_0716/DEEP2_Field'+str(ii)+'_all_line_fit.fits'
    data  = Table(fits.getdata(file1))
    if ii == 1:
        data0 = data
    else:
        data0 = vstack([data0, data])

print 'data0 : ', len(data0)

#Start Graph
fig1, ax1 = plt.subplots() #plt.gcf()

xlim = [0.4,50]
ylim = [0.1,20]

# Gridding of data
R23_bin  = 0.25
O32_bin  = 0.25
R23_grid = np.arange(np.log10(xlim[0]), np.log10(xlim[1])+R23_bin, R23_bin)
O32_grid = np.arange(np.log10(ylim[0]), np.log10(ylim[1])+R23_bin, O32_bin)

#print len(R23_grid), len(O32_grid)

N_arr0 = np.zeros( (len(R23_grid),len(O32_grid)), dtype=np.float)

N = len(R23_grid)
M = len(O32_grid)
print N, M

T_arr = np.zeros((N, M), dtype = object)
print len(T_arr[:,0])

#For Loop that imports data 

O2 = data0['OII_FLUX_MOD']
O3 = 1.33*data0['OIIIR_FLUX_MOD']
Hb = data0['HB_FLUX_MOD']
R23 = (O2+O3)/Hb
O32 = O3/O2
SNR2 = data0['OII_SNR']
SNR3 = data0['OIIIR_SNR']
SNRH = data0['HB_SNR']

#SNR code: This rules out major outliers by only using specified data
det3 = np.where((SNR2 >= 3) & (SNR3 >= 3) & (SNRH >= 3) &
                (O2 > 0) & (O3 > 0) & (Hb>0))[0]

#Plotting
label0 = 'Field'+str(ii)+', N='+str(len(det3))
x = np.log10(R23[det3])
y = np.log10(O32[det3])
finite0 = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
x = x[finite0]
y = y[finite0]
scatter = ax1.scatter(x,y,2.5, facecolor='b', edgecolor='none', marker='o',
                      alpha=0.5, label=label0)
    
    
x0 = x.tolist()
y0 = y.tolist()

for jj in range(len(R23_grid)):
    for kk in range(len(O32_grid)):
        array= np.where((x < R23_grid[jj]+R23_bin) & (x >= R23_grid[jj]) &
                        (y < O32_grid[kk]+O32_bin) & (y >= O32_grid[kk]))[0] 
        N_arr0[jj,kk]    += len(array)
        T_arr[jj,kk] = det3[array]



print np.max(N_arr0), np.min(N_arr0), np.average(N_arr0)

ax1.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for Data Set 1')
ax1.set_xlabel(r'log($R_{23}$)')
ax1.set_ylabel(r'log($O_{32}$)')
ax1.set_xlim(np.log10(xlim))
ax1.set_ylim(np.log10(ylim))
ax1.minorticks_on()
ax1.legend(loc='upper left', numpoints=3) #"Field1", "Field2", "Field3", "Field4")    
fig1.set_size_inches(9,9)
#fig.savefig(fitspath+'R23_O32.pdf')

pdf_pages.savefig(fig1)

fig2 = plt.figure()
ax2 = plt.gca() # ax = plt.subplots() #plt.gcf()
cm= plt.cm.get_cmap('Blues')

#Colorbar and hexbin plotting
tabmastergrid = Table([x0,y0])
asc.write(tabmastergrid, fitspath+'testmastergrid.tbl', format='fixed_width_two_line')
hb0 = ax2.hexbin(x0, y0, gridsize=(len(R23_grid),len(O32_grid)), cmap= cm)
cbaxes= fig2.add_axes([0.135,0.20, 0.75, 0.025])
cb= fig2.colorbar(hb0, cax=cbaxes, ticks=[0.,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300], orientation= 'horizontal')
cb.set_label('density')
#fig.tight_layout()

ax2.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for Data Set 1')
ax2.set_xlabel(r'log($R_{23}$)')
ax2.set_ylabel(r'log($O_{32}$)')
ax2.set_xlim(np.log10(xlim))
ax2.set_ylim(np.log10(ylim))
ax2.minorticks_on()
ax2.legend(loc='upper left', numpoints=3)     

fig2.set_size_inches(8,8)
pdf_pages.savefig(fig2)

pdf_pages.close()

fig1.clear()
fig2.clear()

outfile = 'Arrays_R23O32bin01MasterGrid.npz' 
np.savez(outfile, T_arr=T_arr, R23_grid=R23_grid, O32_grid=O32_grid,
         N_arr0=N_arr0)


