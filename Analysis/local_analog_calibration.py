"""
local_analog_calibration
====

Plots R23 and O32 line ratios and compare it to Bian et al. (2018)
calibration that is based on local analogs to z~2 galaxies
"""

import sys, os

from chun_codes import systime, match_nosort_str

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

def bian18_R23(OH):
    '''
    Function to return log(R23) given metallicity

    Parameters
    ----------
     OH : 12+log(O/H)
    '''

    R23_coeff = [138.0430, -54.8284, 7.2954, -0.32293]
    R23_p = np.poly1d(R23_coeff)

    return R23_p(OH)
#enddef

