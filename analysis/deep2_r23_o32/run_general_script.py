#!/usr/bin/env python3
from Zcalbase_gal.analysis.deep2_r23_o32 import general
from Metallicity_Stack_Commons import get_user

params_list = [
    {'raw': True,  # No MC, no dust, use default valid_table
     'apply_dust': False,
     'revised': False},
    {'raw': True,  # No MC, no dust, use revised valid_table
     'apply_dust': False,
     'revised': True},
    {'raw': True,  # No MC, apply dust, use default valid_table
     'apply_dust': True,
     'revised': False},
    {'raw': True,  # No MC, apply dust, use revised valid_table
     'apply_dust': True,
     'revised': True},
    {'raw': False,  # MC, no dust, use default valid_table
     'apply_dust': False,
     'revised': False},
    {'raw': False,  # MC, no dust, use revised valid_table
     'apply_dust': False,
     'revised': True},
    {'raw': False,  # MC, apply dust, use default valid_table
     'apply_dust': True,
     'revised': False},
    {'raw': False,  # MC, apply dust, use revised valid_table
     'apply_dust': True,
     'revised': True},
]
dataset = 'n_Bins'
general.run_grid_r23_o32_analysis(dataset)

fitspath_ini = get_user()
for params in params_list:
    general.run_grid_plots(fitspath_ini, dataset, **params, individual=False)
