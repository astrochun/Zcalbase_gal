import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages


#ascii_file = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0514_old/bin_info.tbl'
#pdf_file = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0514_old/manual_binning_shading.pdf'

n_split = 3

def graph_bins(ascii_file,pdf_file):
    '''
    Purpose: This function takes inputed grid parameters and graphs them on the log(R23)--log(O32) plane 
    Inputs: 
           ascii_file -> file with the minimum, maximum, and average R23 and O32 value
           pdf_file   -> name of pdf file made by code
    Outputs: 
           none

    Currently works for manual binning with an n_split (see general function) of 3. More updates would need to be made for a different n_split to correct indexing through arrays. 
    '''
    pdf_pages = PdfPages(pdf_file)
    tab = asc.read(ascii_file)

    lR23_min = tab['logR23_min'].data
    lO32_min = tab['logO32_min'].data
    lR23_avg = tab['logR23_avg'].data
    lO32_avg = tab['logO32_avg'].data
    lR23_max = tab['logR23_max'].data
    lO32_max = tab['logO32_max'].data
    #print(lR23_min, type(lR23_min))
    
    
    fig, ax = plt.subplots()
    #finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    #x1 = R23[finite0]
    #y1 = O32[finite0]
    #x = np.log10(x1)
    #y = np.log10(y1)
    #vlines = np.log10(R23_lowlimit)
    #hlines = np.log10(O32_lowlimit)
    #ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for Manual Binning')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    #for pp in range(n_split*len(bin_start)): plt.axvline(x = vlines[pp], linewidth= 0.3, color= 'k')

    
    count = 0
    for aa in range(len(lR23_min)):
        xmin = lR23_min[aa]
        if aa <= (len(lO32_min)-n_split-1):
            xmax = lR23_max[aa]
        else:
            xmax = lR23_max[-1]
        #plt.axvline(x = lR23_min[jj], linewidth= 0.3, color= 'k')
        
        x_value = [xmin,xmax]
        y_value = [lO32_min[aa], lO32_min[aa]]
        y_average = [lO32_avg[aa],lO32_avg[aa]]
        #plt.plot(x_value,y_value, linewidth= 0.3, color= 'b')
        plt.plot(x_value, y_average, linewidth= 0.75, color= 'k')


        #for cc in range(3):
        if aa<=(len(lR23_min)-2):
            xx= [xmin, lR23_max[aa]]
        else:
            xx = [lR23_min[aa-1],lR23_max[-1]]
            #print('xx= ', xx)  
        if aa <=(len(lO32_max)-1):
            ymax = lO32_max[aa]
        else:
            ymax = lO32_max[-1]
        print('ymin, ymax = ', lO32_min[aa], ymax)
        ax.fill_between(xx, lO32_min[aa], ymax)
        count+=1
    print(count)
    fig.savefig(pdf_pages, format ='pdf')
    pdf_pages.close()

    
        
       
    
