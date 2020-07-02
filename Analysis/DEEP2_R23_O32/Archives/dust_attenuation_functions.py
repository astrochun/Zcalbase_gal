
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


'''k_ascii = fitspath_ini+'/cardelli_k_values.tbl'
    n3 = ('Line_Name','K_value')
    tab3 = Table([line_name, k_values])
    asc.write(tab3, k_ascii, format = 'fixed_width_two_line')'''




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

