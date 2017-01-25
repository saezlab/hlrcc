#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, Series
from framed import load_cbmodel, simplify, FBA, MOMA, pFBA, FVA, lMOMA
from hlrcc.metabolomics.sampler import sample

# -- Imports
# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.csv', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .05]

# Medium
medium = read_csv('./files/DMEM_41966_medium.txt', sep='\t').dropna().set_index('exchange')

# O2 consumption (mmol/gDW/h)
o2_exch = 'R_EX_o2_e'
core_o2 = read_csv('./data/uok262_metabolomics_core_o2_processed.txt', sep='\t', index_col=0)[o2_exch]
core_o2 *= -1

# -- Import metabolic model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# FVA to simplify model
simplify(model)
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# -- Fit medium
conditions = ['UOK262', 'UOK262pFH']

res_fba, res_env = {}, {}

# c = 'UOK262'
for c in conditions:
    print c
    res_fba[c] = {}

    # MOMA with environmental conditions restricted to metabolites in the medium
    env = {r: (-10, model.reactions[r].ub) if r in list(medium.index) else (0, model.reactions[r].ub) for r in model.reactions if r.startswith('R_EX_') or r.startswith('R_sink_')}

    # Add FH KO to tumour cell lines
    if c == 'UOK262':
        env.update({'R_FUM': (0, 0), 'R_FUMm': (0, 0)})

    # Minimise differences of measured rates: [(r, meas[r], moma_sol.values[r]) for r in meas]
    meas = {r: core.ix[r, c] if r != o2_exch else core_o2[c] for r in list(core.index) + [o2_exch] if r in model.reactions}

    moma_sol = lMOMA(model, reference=meas, constraints=env)
    env.update({r: moma_sol.values[r] for r in meas})
    print moma_sol.fobj

    # # [(r, env[r]) for r in env if isinstance(env[r], float) or env[r][0] != 0]
    # meas_medium = {r: 0 for r in medium.index if r in model.reactions and r not in meas}
    # moma_medium_sol = lMOMA(model, reference=meas_medium, constraints=env)
    # env.update({r: moma_medium_sol.values[r] for r in meas_medium})
    # print moma_medium_sol.fobj

    res_env[c] = env

    # Biomass and ATP production
    res_fba[c]['biomass'] = FBA(model, objective={'R_biomass_reaction': 1}, constraints=env).fobj
    print res_fba[c]['biomass']

    res_fba[c]['atp'] = FBA(model, objective={'R_ATPM': 1}, constraints=env, get_shadow_prices=True).fobj
    print res_fba[c]['atp']

    # Sampler
    sampling = sample(model, n_samples=100, n_steps=200, verbose=1, constraints=env)
    sampling.to_csv('./data/%s_sampling.txt' % c, sep='\t', index=False)
    print '[INFO] Sampling finished: ', c

# res_fba = DataFrame(res_fba).T
#
# print 'ATP yield / mol of Glucose: UOK262pFH %.2f, UOK262 %.2f' % (
#     (res_fba.ix['UOK262pFH', 'atp'] / -res_env['UOK262pFH']['R_EX_glc_e']), (res_fba.ix['UOK262', 'atp'] / -res_env['UOK262']['R_EX_glc_e'])
# )
#
# print 'Biomass yield / mol of Glucose: UOK262pFH %.2f, UOK262 %.2f' % (
#     (res_fba.ix['UOK262pFH', 'biomass'] / -res_env['UOK262pFH']['R_EX_glc_e']), (res_fba.ix['UOK262', 'biomass'] / -res_env['UOK262']['R_EX_glc_e'])
# )
