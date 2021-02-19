import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii as asc

fitspath_ini = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'
temp_met = '/Users/reagenleimbach/Desktop/Zcalbase_gal/R23O32_Manual_0202/n_Bins_temperatures_metalicity.tbl'  #fitspath_ini+'R23O32_Manual_0202/n_Bins_temperature_metalicity.tbl'
tm_table = asc.read(temp_met)

valid_file = fitspath_ini + 'verification_tables/n_Bins_verification_tab.tbl'
valid_tab = asc.read(valid_file)

out_pdf = fitspath_ini+'green_pea_jiang.pdf'

# This is a keyword that will be passed in in order to specify which calibration we want to compare to
y_calibration = 'jiang'

O32_all = tm_table['O32_Composite']
R23_all = tm_table['R23_Composite']
com_O_log = tm_table['com_O_log']  # This is the 12+log(OH) value
ID = tm_table['ID']

detect = valid_tab['Detection']

det_4363 = np.where((detect == 1.0))[0]
nan_det_4363 = np.where((detect == 0.0))[0]
r_limit = np.where((detect == 0.5))[0]

det_O32 = O32_all[det_4363]
det_R23 = R23_all[det_4363]
det_OH = com_O_log[det_4363]
det_ID = ID[det_4363]
    
nandet_O32 = O32_all[nan_det_4363]
nandet_R23 = R23_all[nan_det_4363]
print(len(nandet_O32), len(nandet_R23))
nandet_OH = com_O_log[nan_det_4363]
nandet_ID = ID[nan_det_4363]


lR23 = det_R23
lO32 = det_O32
OH = det_OH

# green_peas_calibration.main(lR23, lO32, OH, out_pdf, n_bins = 6, xra = [0.3,1.5], yra = [6.5,9.10],
# marker=['D','X','3','4'], label=['Detection','Non-Dectection','DEEP2', 'MACT'], fit=False, silent=False, verbose=True)

# Given R23_composite and O32_composite, calculate the metallicity according to the Jiang calibration,
# and plot against Zcalbase_gal derived metallicities
a = -24.135
b = 6.1523
c = -0.37866
d = -0.147
e = -7.071

x_plus = np.zeros(len(det_R23))
x_minus = np.zeros(len(det_R23))
gamma = np.zeros(len(det_R23))

for ii in range(len(det_R23)): 
    # gamma[ii] = -(lR23[ii]- a+d*e*lO32[ii])
    gamma[ii] = np.abs((b-d*lO32[ii])*(b-d*lO32[ii])-4*c*(-(lR23[ii]- a+d*e*lO32[ii])))
    x_plus[ii] = (-(b-d*lO32[ii])+np.sqrt(gamma[ii]))/(2*c)
    x_minus[ii] = (-(b-d*lO32[ii])-np.sqrt(gamma[ii]))/(2*c)

print(det_ID)

fig, ax = plt.subplots()
ax.scatter(x_plus, det_OH, c='b')
ax.scatter(x_minus, det_OH, c='g')
ax.set_xlabel('Jiang Calculation')
ax.set_ylabel('Composite Metallicities')

for qq in range(len(det_ID)):
    ax.annotate(det_ID[qq], (x_plus[qq], det_OH[qq]), fontsize='12')
for aa in range(len(det_ID)):
    ax.annotate(det_ID[aa], (x_minus[aa],det_OH[aa]), fontsize='12')
