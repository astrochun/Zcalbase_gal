#!/usr/bin/env python3
from Zcalbase_gal.analysis.deep2_r23_o32 import general
from Metallicity_Stack_Commons import get_user
dataset = 'n_Bins'
general.run_grid_r23_o32_analysis(dataset)

fitspath_ini = get_user()
general.run_grid_plots(fitspath_ini, dataset, raw=False, apply_dust=False,
                       revised=False, individual=False)
general.run_grid_plots(fitspath_ini, dataset, raw=False, apply_dust=False,
                       revised=True, individual=False)
general.run_grid_plots(fitspath_ini, dataset, raw=True, apply_dust=False,
                       revised=False, individual=False)
general.run_grid_plots(fitspath_ini, dataset, raw=True, apply_dust=False,
                       revised=True, individual=False)
