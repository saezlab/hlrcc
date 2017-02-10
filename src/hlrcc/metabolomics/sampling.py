#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
from pandas import read_csv, DataFrame
from hlrcc.metabolomics.sampler import sample, fix_futile_cycles
from framed import load_cbmodel, simplify, pFBA, lMOMA


# -- Import metabolic model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))


# -- Imports
conditions = ['UOK262', 'UOK262pFH']

# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.csv', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .1]

# Medium
medium = read_csv('./files/DMEM_41966_medium_revised.txt', sep='\t').dropna().set_index('exchange')

# Doubling-times and biomass calculation
dt = read_csv('./data/core/uok262_doubling_times.csv')
dt = dt[dt['experiment'] != 'rep3']
growthrate = {c: np.float(np.mean(np.log(2) / dt.loc[dt['condition'] == c, 'DW doubling time (h)'])) for c in conditions}


# -- Fit medium
res_fba = {}

# c = 'UOK262'
for c in conditions:
    print c
    res_fba[c] = {}
    c_model = model.copy()

    # MOMA with environmental conditions restricted to metabolites in the medium
    env = {r: (medium.ix[r, 'lb'], c_model.reactions[r].ub) if r in list(medium.index) else (0, c_model.reactions[r].ub) for r in c_model.reactions if r.startswith('R_EX_') or r.startswith('R_sink_') or r.startswith('R_DM_')}

    # Add FH KO to tumour cell lines
    if c == 'UOK262':
        env.update({'R_FUM': (0, 0), 'R_FUMm': (0, 0)})

    # Fix Biomass
    env.update({'R_biomass_reaction': growthrate[c]})

    # Minimise differences of measured rates: [(r, meas[r], moma_sol.values[r]) for r in meas]
    meas = {r: core.ix[r, c] for r in list(core.index) if r in c_model.reactions}
    moma_sol = lMOMA(c_model, reference=meas, constraints=env)
    env.update({r: moma_sol.values[r] for r in meas})
    print moma_sol

    # Write environmental conditions in the model
    for rid in env:
        if type(env[rid]) is float:
            c_model.set_flux_bounds(rid, env[rid], env[rid])
        else:
            c_model.set_flux_bounds(rid, env[rid][0], env[rid][1])

    # FVA to simplify model
    simplify(c_model)
    print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(c_model.metabolites), len(c_model.reactions), len(c_model.genes))

    # Biomass and ATP production
    res_fba[c]['biomass'] = pFBA(c_model, objective={'R_biomass_reaction': 1})
    print 'Max biomass: ', res_fba[c]['biomass'].pre_solution.fobj

    res_fba[c]['atp'] = pFBA(c_model, objective={'R_ATPM': 1})
    print 'Max ATP: ', res_fba[c]['atp'].pre_solution.fobj

    # Sampling
    for sample_rep in range(1, 3):
        samples = sample(c_model, n_samples=1000, n_steps=2500, verbose=1)
        samples.to_csv('./data/%s_sampling_%d.csv' % (c, sample_rep), index=False)

        samples = fix_futile_cycles(c_model, samples=samples)
        samples.to_csv('./data/%s_sampling_processed_%d.csv' % (c, sample_rep), index=False)
        print '[INFO] Sampling finished: ', c

# -- Store fluxes
fluxes = DataFrame({c: res_fba[c]['atp'].values for c in res_fba})
fluxes['delta'] = fluxes['UOK262'].abs() - fluxes['UOK262pFH'].abs()
fluxes.to_csv('./data/pfba_atp.csv')
