
###THIS BINS DATA FOR ALL BINNING METHODS AND IS USED IN GENERAL FUNCTION###



#Numpy binary file or FITS binary table: read up on both and attempt both
#Change name of pdf, name of output file, bin size

###Create r 23 grid to have nxm dimensions and update each of the corresponding values, like a table
### m will be 2 and



import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table


def n_times_binned(fitspath, pdf_pages, outfile, n_split, individual_ID, R23,O32,
                   O2, O3, Hb, SNR2, SNR3, SNRH, det3, data3,galinbin, adaptive=False):
    dataset = 'n_Bins'
    ##R23 and O32 are going to be log values --> check to make sure that this doesn't affect the binning
    #One_dimensional binning for R23 followed by each bin being split in O32 in n_split bins 
    '''increase the number of galaxies as R23 increases'''
    pdf_pages = PdfPages(pdf_pages)

    #One_dimensional binning for R23 

    #R23 = 10**(logR23)
    #O32 = 10**(logO32)
    sortR23 = np.argsort(R23)
    R23_sort0 = R23[sortR23]


    n_bins = len(galinbin)
    
    #print(n_bins, 'n_split',n_split)
    n_bins_range = np.arange(0,n_split*n_bins,1)


    #Initializing Arrays for Grid stacking
    Number_inbin = np.zeros((n_bins,n_split))   #Used to be N_arr0
    locator  = np.zeros((n_bins,n_split), dtype = object) #Used to be T_arr
    O32_minimum    = np.zeros((n_bins,n_split))    #Used to be O32_grid
    R23_minimum    = np.zeros((n_bins,n_split))      #Used to be R23_grid
    O32_median = np.zeros((n_bins,n_split))
    R23_median = np.zeros((n_bins,n_split))
    O32_max    = np.zeros((n_bins,n_split))
    R23_max    = np.zeros((n_bins,n_split))
    xBar  = np.zeros(n_split*n_bins)
    yBar  = np.zeros(n_split*n_bins)
    area  = np.zeros(n_split*n_bins)
    Bin_number = np.zeros(len(data3), dtype = int) #Used to be N_bin

    R23_minall = np.zeros(len(R23))
    O32_minall = np.zeros(len(R23))
    R23_avgall = np.zeros(len(R23))
    O32_avgall = np.zeros(len(R23))
    R23_medall = np.zeros(len(R23))
    O32_medall = np.zeros(len(R23))
    R23_maxall = np.zeros(len(R23))
    O32_maxall = np.zeros(len(R23))

    #Bin starts and stops initializing
    bin_start = np.zeros(n_bins, dtype=np.int)
    bin_end   = np.zeros(n_bins, dtype=np.int)

    for ii in range(n_bins):
        print('ii: ', ii) 
        if ii == 0:
            bin_start[ii] = 0
        else:
            bin_start[ii]= bin_end[ii-1]+1
        if ii == n_bins-1:
            bin_end[ii] = len(R23_sort0)-1
        else:
            bin_end[ii] = galinbin[ii]+bin_start[ii]-1   
        print('Bin Start:', bin_start[ii] , 'Bin end:', bin_end[ii])


        
        R23_idx = np.where((R23>= R23_sort0[bin_start[ii]]) & (R23<= R23_sort0[bin_end[ii]]))[0]
        # There is probably a more streamline way to do this but let's start with version 1 first then update later
        # @{

        R23_inbins = R23[R23_idx]
        O32_inbins = O32[R23_idx]          #O32 relative to the R23_index/ O32 sorted into the R23 bins
        O32_index = np.argsort(O32_inbins) #Sort the O32 in their bins so that we have their locations
        sortO32 = O32_inbins[O32_index]    #Take the O32 in each bin and organizes them based on the index in the previous line

        # }

        # The following lines could go into a defintion called "evenSplit( ... )" therefore
        # allowing you to create a different splitting method definition called "optimalSplit( ... )"
        # or something to that effect if so desired.
        # @{
        n_subbins = np.int(np.floor(float(len(sortO32))/n_split))
        subbin_arr = np.ones(n_split, dtype = int)*n_subbins #We have to take the length of the range of n_split because n_split is of type int
        n_remainder = len(sortO32) - n_subbins*n_split
        backwardsIdx = -1
        forwardsIdx = 0
        
        # Sorting all the remainders into the correct bins: Use backwardsIdx or forwardsIdx to start adding spectra
        #from the front or the back of the list of bins 
        for rr in range(n_remainder):
            subbin_arr[backwardsIdx] = n_subbins+1
            backwardsIdx -= 1
        subbin_arr[backwardsIdx] = n_subbins
        # }
        
        startIdx = 0
        endIdx = subbin_arr[0]
        for jj in range(n_split):
            # Let's grab all O32 values
            O32_values_perbin = sortO32[startIdx:endIdx]
            O32_inbins_idx = O32_index[startIdx:endIdx]

            #This index gives the positions of all the R23 values relative 
            #to the O32 values over the entire R23 bin
            N_bin_idx = R23_idx[O32_inbins_idx]        

            #Gives the bin number for each spectra
            Bin_number[N_bin_idx] = (ii*n_split)+jj  

            # Now let's start storting our data into variables to use later

            #First two map minimum R23 and O32 measure to all individual spectra in bin 
            #Second two values give the lowest O32 and R23 value set for each bin
            R23_minall[N_bin_idx]= R23_sort0[bin_start[ii]]
            O32_minall[N_bin_idx]= O32_values_perbin[0]
            R23_minimum[ii,jj]   = R23_sort0[bin_start[ii]]
            O32_minimum[ii,jj]   = O32_values_perbin[0]

            #First two map median R23 and O32 measure to all individual spectra in bin 
            #Second two give the median R23 and O32 value for each bin
            R23_medall[N_bin_idx]= np.median(R23[N_bin_idx])
            O32_medall[N_bin_idx]= np.median(O32_values_perbin)      
            R23_median[ii,jj]    = np.median(R23[N_bin_idx])
            O32_median[ii,jj]    = np.median(O32_values_perbin)     

            #First two map maximum R23 and O32 measure to all individual spectra in bin 
            #Second two give the maximum R23 and O32 value for each bin
            R23_maxall[N_bin_idx]= np.max(R23[N_bin_idx])
            O32_maxall[N_bin_idx]= np.max(O32[N_bin_idx])
            R23_max[ii,jj]       = np.max(R23[N_bin_idx])
            O32_max[ii,jj]       = np.max(O32[N_bin_idx])

            #Gives the number of galaxies in each bin
            Number_inbin[ii,jj]  += len(O32_values_perbin)                      

            #Gives the index (location numbers) for the spectra in each bin
            #and is used later loop over and get the galaxies 
            locator[ii,jj]   = N_bin_idx                                    

            #Maps average R23 measure to all individual spectra in bin 
            R23_avgall[N_bin_idx]= np.average(R23[N_bin_idx])
            #Gives the average R23 value for each bin
            xBar[(ii*n_split)+jj]= np.average(R23[N_bin_idx])

            #Maps average O32 measure to all individual spectra in bin
            O32_avgall[N_bin_idx]= np.average(O32_values_perbin)
            #Gives the average O32 value for each bin
            yBar[(ii*n_split)+jj]= np.average(O32_values_perbin)

            #Gives the number of galaxies in each bin
            area[(ii*n_split)+jj]= len(O32_values_perbin)                 

            # Now we can shift our window over for the next split like the counting index the Stacking code
            startIdx = endIdx
            endIdx = startIdx+subbin_arr[jj]+1

    O32_lowlimit = O32_minimum.reshape(n_split*n_bins) #Used to be O32_grids
    R23_lowlimit = R23_minimum.reshape(n_split*n_bins) #Used to be R23_grids
    O32_medians = O32_median.reshape(n_split*n_bins)
    R23_medians = R23_median.reshape(n_split*n_bins)
    R23_maxval = R23_max.reshape(n_split*n_bins)
    O32_maxval = O32_max.reshape(n_split*n_bins)

    #Plotting
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    vlines = np.log10(R23_lowlimit)
    hlines = np.log10(O32_lowlimit)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')
    #for pp in range(n_split*len(bin_start)): plt.axvline(x = vlines[pp], linewidth= 0.3, color= 'k')

    for jj in range(len(O32_lowlimit)):
        xmin = vlines[jj]
        if jj <= (len(O32_lowlimit)-n_split-1): xmax = vlines[jj+n_split]
        else: xmax = np.log10(max(R23))
        plt.axvline(x = vlines[jj], linewidth= 0.3, color= 'k')
        
        x_value = [xmin,xmax]
        y_value = [hlines[jj], hlines[jj]]
        y_average = [yBar[jj],yBar[jj]]
        plt.plot(x_value,y_value, linewidth= 0.3, color= 'b')
        plt.plot(x_value, y_average, linewidth= 0.3, color= 'g')


        #plt.xlim(-0.3, 1)
    fig.savefig(pdf_pages, format ='pdf')


    '''
    fig, ax = plt.subplots()
    finite0 = np.where((np.isfinite(R23)) & (np.isfinite(O32)))[0]
    x1 = R23[finite0]
    y1 = O32[finite0]
    x = np.log10(x1)
    y = np.log10(y1)
    vlines = np.log10(R23_lowlimit)
    hlines = np.log10(O32_lowlimit)
    ax.scatter(x,y,1.5, facecolor='r', edgecolor='face', marker='*',alpha=1)
    ax.set_title(r'$R_{23}$ vs. $O_{32}$ Plot for DEEP2')
    ax.set_xlabel(r'log($R_{23}$)')
    ax.set_ylabel(r'log($O_{32}$)')

    for aa in range(len(vlines)-1):
        if aa < (len(vlines)-1):
            x = np.linspace(vlines[aa],vlines[aa+1])
        else: x = np.linespace(
    ax.fill_between()

    for jj in range(len(O32_lowlimit)):
        xmin = vlines[jj]
        if jj <= (len(O32_lowlimit)-n_split-1): xmax = vlines[jj+n_split]
        else: xmax = np.log10(max(R23))
        #print "jj, xmin, xmax, hlines[jj], vlines[jj]", jj, xmin, xmax, hlines[jj], vlines[jj]
        plt.axvline(x = vlines[jj], linewidth= 0.3, color= 'k')
        
        x_value = [xmin,xmax]
        y_value = [hlines[jj], hlines[jj]]
        y_average = [yBar[jj],yBar[jj]]
        plt.plot(x_value,y_value, linewidth= 0.3, color= 'b')
        plt.plot(x_value, y_average, linewidth= 0.3, color= 'g')



        #plt.xlim(-0.3, 1)
    fig.savefig(pdf_pages, format ='pdf')'''

    pdf_pages.close()


    np.savez(outfile, locator=locator, R23_minimum=R23_minimum, O32_minimum=O32_minimum, Number_inbin=Number_inbin)
    
    n1 = ('bin_ID' ,'N_stack','logR23_min','logO32_min','logR23_avg','logO32_avg',
          'logR23_median','logO32_median','logR23_max','logO32_max')
    tab1 = Table([n_bins_range,area,np.log10(R23_lowlimit), np.log10(O32_lowlimit),np.log10(xBar),
                  np.log10(yBar),np.log10(R23_medians),np.log10(O32_medians),np.log10(R23_maxval),
                  np.log10(O32_maxval)], names = n1)
    asc.write(tab1, fitspath+'/bin_info.tbl', format='fixed_width_two_line')   #used to be called +dataset+'_binning_averages.tbl
    

    '''variable_formats = {'bin_ID':' %i ', 'N_stack': '%i', 'logR23_min': '%.2f', 'logO32_min':  '%.2f', 'logR23_avg':  '%.2f',
                        'logO32_avg':  '%.2f','logR23_median':  '%.2f','logO32_median':  '%.2f','logR23_max':  '%.2f','logO32_max':  '%.2f'}
    asc.write(tab1, '/Users/reagenleimbach/Desktop/Zcalbase_gal/Honors_Thesis/mostrecent_bin_info.tex', format='latex', formats= variable_formats)
    
    fig.clear()'''

    n2=('logR23', 'logO32', 'OIII_5007_S/N', 'bin_ID','ID',
        'logR23_min', 'logO32_min', 'logR23_avg','logO32_avg',
        'logR23_median','logO32_median','logR23_max','logO32_max')
    tab2= Table([R23, O32, SNR3, Bin_number, individual_ID,
                 R23_minall,O32_minall,R23_avgall,O32_avgall,
                 R23_medall,O32_medall,R23_maxall, O32_maxall], names=n2)
    asc.write(tab2, fitspath+'/individual_bin_info.tbl', format='fixed_width_two_line') #used to be + dataset+'_2d_binning_datadet3.tbl

    '''n3 = ('ID' , 'R23_grid', 'O32_grid')
    tab1 = Table([n_bins_range, R23_grid, O32_grid], names = n3)
    asc.write(tab1, fitspath+'/Double_Bin_grid_values.tbl', format='fixed_width_two_line')'''
###Create another ascii table with the R23_grid and O32_grid values for plots 
