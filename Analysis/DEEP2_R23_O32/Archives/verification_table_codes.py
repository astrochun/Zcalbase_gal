import numpy as np
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import vstack, hstack
from astropy.table import Table
import os
from os.path import exists
import glob
from datetime import date



def check_verification_table(fitspath_ini, dataset, combine_flux_ascii):
    verification_table = fitspath_ini+'verification_tables/'+dataset+'_verification_tab.tbl'
    if exists(verification_table):
        return verification_table
    else:
        print('Making verification table')
        verification_tables.verification_master(fitspath,dataset, combine_flux_ascii)
        return verification_table
