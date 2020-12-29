
def dust_attenuation(fitspath, combine_ascii):
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    
    combine_asc = asc.read(combine_ascii)
    ini_con = 0.468
    ID = combine_asc['ID']
    HBeta= combine_asc['HBETA_Flux_Observed'].data
    HGamma= combine_asc['HGAMMA_Flux_Observed'].data
    
    lam0_OII = combine_asc['OII_3727_X_bar'].data
    lam0_HDELTA = combine_asc['HDELTA_X_bar'].data
    lam0_Hgamma = combine_asc['HGAMMA_X_bar'].data
    lam0_HBETA = combine_asc['HBETA_X_bar'].data
    lam0_4363 = combine_asc['OIII_4363_X_bar'].data
    lam0_4958 = combine_asc['OIII_4958_X_bar'].data
    lam0_5007 = combine_asc['OIII_5007_X_bar'].data

    k_3727 = call_cardelli(lam0_OII)
    k_HDELTA = call_cardelli(lam0_HDELTA)
    k_Hgamma = call_cardelli(lam0_Hgamma)
    k_HBETA = call_cardelli(lam0_HBETA)
    k_4363 = call_cardelli(lam0_4363)
    k_4958 = call_cardelli(lam0_4958)
    k_5007 = call_cardelli(lam0_5007)

    
    
    EBV= np.log10((HBeta/HGamma)*(ini_con))*2.5*(1/(k_Hgamma-k_HBETA))  
    for nn in range(len(HGamma)):
        if EBV[nn] <= 0: EBV[nn] = 0
    
    #print EBV
    A_3727 = EBV*k_3727
    A_HDELTA = EBV*k_HDELTA
    A_Hgamma = EBV*k_Hgamma
    A_HBETA = EBV*k_HBETA
    A_4363 = EBV*k_4363
    A_4958 = EBV*k_4958
    A_5007 = EBV*k_5007
    #print "A_3727:", A_3727

    out_ascii = fitspath+'/dust_attentuation_values.tbl'
    #if not exists(out_ascii_single):
    n2= ('ID','k_3727', 'k_HDELTA', 'k_Hgamma', 'k_HBETA', 'k_4363', 'k_4958', 'k_5007', 'E(B-V)')
    tab1 = Table([ID,k_3727, k_HDELTA, k_Hgamma, k_HBETA , k_4363, k_4958, k_5007, EBV], names=n2)
    asc.write(tab1, out_ascii, format='fixed_width_two_line')
    
    
def call_cardelli(lam0): #, extrapolate=False):
    #lambda0 =[3726.16, 3868.74, 3888.65, 3967.51, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]* u.angstrom
    line_name = ['OII_3727','NeIII','HeI','3967', 'HDELTA', 'Hgamma', 'OIII_4363', 'HBETA', 'OIII_4958','OIII_5007']
    lambda0 = lam0*u.angstrom
    k_values= cardelli(lambda0,R=3.1)
    return k_values




def dust_vs_nondust_table(fitspath, dust_metal_table, nondust_metal_table, dust_atten_values, name):
    #dust_vs_nondust_table = fitspath + 'dust_and_nondust_metaltab.tbl'
    dust_vs_nondust_table = fitspath + name
    
    #Non Dust attentuation
    nondust = asc.read(nondust_metal_table)
    dust = asc.read(dust_atten_values)
    
    nondust_metal = nondust['com_O_log'].data
    dust_metal = dust['com_O_log'].data
    ID = nondust['ID'].data
    R23_composite = nondust['R23_Composite'].data
    O32_composite = nondust['O32_Composite'].data
    Temperature = nondust['Temperature'].data

    n_dust = ('ID','R23_Composite','O32_Composite', 'Non-Dust Attenuated Metallicities','Dust Attenuated Metallicities','Temperature')
    tab_dust = Table([ID, R23_composite, O32_composite, nondust_metal, dust_metal, Temperature], names =n_dust)
    asc.write(tab_dust, dust_vs_nondust_table, format = 'fixed_width_two_line')


def dust_att_plot(fitspath, combine_flux):
    """
    Purpose
    Produces pdf file of plots comparing the average of HGAMMA/HBETA of each bin to R23 and O32
    Preliminary plots for dust attenuation work
    """
    pdf_pages = PdfPages(join(fitspath, 'dust_attenuation_plots.pdf'))
    com_asc = asc.read(combine_flux)
    H_gamma_obs = com_asc['Hgamma_Flux_Observed']
    H_beta_obs = com_asc['OIII_4958_Flux_Observed']
    R23 = com_asc['R_23_Average']
    O32 = com_asc['O_32_Average']

    Gambet = H_gamma_obs / H_beta_obs

    fig, ax = plt.subplots()
    ax.scatter(Gambet, O32, marker='.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('O32')
    ax.set_title('H_gamma/H_beta vs. O32')
    fig.set_size_inches(8, 8)
    fig.savefig(pdf_pages, format='pdf')

    fig, ax = plt.subplots()
    ax.scatter(Gambet, R23, marker='.')
    ax.set_xlabel('H_gamma/H_beta')
    ax.set_ylabel('R23')
    ax.set_title('H_gamma/H_beta vs. R23')
    fig.set_size_inches(8, 8)
    fig.savefig(pdf_pages, format='pdf')
    pdf_pages.close()
