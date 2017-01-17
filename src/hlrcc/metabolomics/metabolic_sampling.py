#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import time
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from pandas import read_csv, DataFrame, Series
from pymist.simulation import min_differences
from hlrcc.metabolomics.sampler import sample
from framed import load_cbmodel, simplify, FVA, FBA

# -- Imports
# Metabolic model
model = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml')
model.remove_b_metabolites()

# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.txt', sep='\t', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .05]

# Medium
medium = read_csv('./files/DMEM_41966_medium.txt', sep='\t').dropna()

# O2 consumption (umol/ugDW/h)
o2_exch = 'R_EX_o2_e_'

core_o2 = read_csv('./data/uok262_metabolomics_core_o2_processed.txt', sep='\t', index_col=0)[o2_exch]
core_o2 /= 1000
core_o2 = core_o2.to_dict()

# -- Fit medium
conditions = ['UOK262', 'UOK262pFH']

for condition in conditions:
    # condition = conditions[0]
    print '[INFO] Condition: ', condition

    # -- Set model constraint variables
    exchange_rates = core[condition].dropna()

    # metabolite = exchange_rates.index[0]
    # Reactions not used
    for reaction in ['R_EX_retpalm', 'R_EX_retpalm_e_']:
        model.set_constraint(reaction, 0)

    # Deepcopy metabolic model
    reduced_model = model.deepcopy()

    # Add FH KO to tumour cell lines
    if condition == 'UOK262':
        reduced_model.set_constraint('R_FUM', 0, 0)
        reduced_model.set_constraint('R_FUMm', 0, 0)
        print '[INFO] R_FUM blocked'

    # -- Fit CORE and medium
    # Block exhange fluxes not present in the medium or CORE experiment
    for reaction in reduced_model.get_exchanges(check_matrix=True):
        if (reaction not in exchange_rates.index) and (reaction != o2_exch) and (reaction not in medium['exchange'].values):
            reduced_model.set_constraint(reaction, lower_bound=0)
    print '[INFO] Exchange reactions constrained'

    # Minimise differences between simulations and measurements
    constrains = {}
    for r in reduced_model.get_exchanges(check_matrix=True):
        if r == o2_exch:
            constrains[r] = -core_o2[condition]
        elif r in exchange_rates.index:
            constrains[r] = exchange_rates.ix[r, condition]

    fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model, check_matrix=True)
    print Series(fitted_medium).sort_values()

    # Constraint metabolic model lower bound
    for reaction in fitted_medium:
        if (reaction in exchange_rates.index) or (reaction == o2_exch):
            reduced_model.set_constraint(reaction, lower_bound=fitted_medium[reaction], upper_bound=fitted_medium[reaction])

    # # -- Metabolic sampling
    # samples = sample(reduced_model, 200, 2500, verbose=1)
    # samples = DataFrame(samples, columns=reduced_model.reactions.keys())
    #
    # # -- Export sampling
    # samples.to_csv('./data/%s_sampling.txt' % condition, sep='\t', index=False)
    # print '[INFO] Sampling finished: ', condition
