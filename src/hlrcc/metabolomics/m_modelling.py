#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv
from framed.experimental.medium import minimal_medium
from framed import load_cbmodel, simplify, FVA, FBA
from framed.cobra.plotting import plot_flux_envelope

medium = read_csv('./files/DMEM_41966_medium.txt', sep='\t').dropna().set_index('exchange')
atp_yileds = read_csv('./files/recon2.2_atp_yield.csv').dropna().set_index('exchange')

# -- Metabolic modelling
# Load
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

s_model = simplify(model, inplace=False)
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(s_model.metabolites), len(s_model.reactions), len(s_model.genes))

pred_atp_yields = {}
for aa in atp_yileds.index:
    env = {r: (0, s_model.reactions[r].ub) for r in s_model.reactions if r.startswith('R_sink_') or r.startswith('R_EX_')}
    env[aa] = (-1, 0)
    env['R_EX_o2_e'] = (None, s_model.reactions['R_EX_o2_e'].ub)

    solution = FBA(s_model, objective={'R_ATPM': 1}, constraints=env)
    print aa, solution.fobj, atp_yileds.ix[aa, 'yield']
    pred_atp_yields[aa] = solution.fobj
