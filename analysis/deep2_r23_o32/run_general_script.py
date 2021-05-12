#!/usr/bin/env python3

from Zcalbase_gal.analysis.deep2_r23_o32 import general
from Metallicity_Stack_Commons import get_user

params_list = [
    # No MC, no dust, use revised valid_table
    {'raw': True, 'apply_dust': False, 'revised': True},
    # No MC, apply dust, use revised valid_table
    {'raw': True, 'apply_dust': True, 'revised': True},
    # MC, no dust, use revised valid_table
    {'raw': False, 'apply_dust': False, 'revised': True},
    # MC, apply dust, use revised valid_table
    {'raw': False, 'apply_dust': True, 'revised': True},
]


if __name__ == '__main__':
    dataset = 'n_Bins'

    general.run_grid_r23_o32_analysis(dataset, apply_dust=True)

    fitspath_ini = get_user()
    for params in params_list:
        general.run_grid_plots(fitspath_ini, dataset, **params,
                               individual=False)
