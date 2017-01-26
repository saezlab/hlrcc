#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, Series
from hlrcc.metabolomics.sampler import sample
from framed import load_cbmodel, simplify, FBA, MOMA, pFBA, FVA, lMOMA, reaction_deletion, gene_deletion

# -- Imports
# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.csv', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .1]

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

# c = 'UOK262pFH'
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
    res_fba[c]['biomass'] = pFBA(model, objective={'R_biomass_reaction': 1}, constraints=env)
    print res_fba[c]['biomass'].pre_solution.fobj

    res_fba[c]['atp'] = pFBA(model, objective={'R_ATPM': 1}, constraints=env)
    print res_fba[c]['atp'].pre_solution.fobj

    # solution = res_fba[c]['atp']
    # print solution.show_metabolite_balance('M_succoa_m', model)
    # print solution.values['R_SUCOASm']

    # # Sampler
    # env['R_ATPM'] = res_fba[c]['atp'].pre_solution.fobj
    #
    # sampling = sample(model, n_samples=200, n_steps=200, verbose=1, constraints=env)
    # sampling.to_csv('./data/%s_sampling.txt' % c, sep='\t', index=False)
    # print '[INFO] Sampling finished: ', c


# -- Plotting
pal = dict(zip(*(conditions, sns.light_palette('#34495e', 3, reverse=True).as_hex()[:-1])))

# ATP and biomass yield
plot_df = DataFrame([
    {
        'condition': 'UOK262pFH',
        'atp': abs(res_fba['UOK262pFH']['atp'].pre_solution.fobj / res_env['UOK262pFH']['R_EX_glc_e']),
        'biomass': abs(res_fba['UOK262pFH']['biomass'].pre_solution.fobj / res_env['UOK262pFH']['R_EX_glc_e'])
    },
    {
        'condition': 'UOK262',
        'atp': abs(res_fba['UOK262']['atp'].pre_solution.fobj / res_env['UOK262']['R_EX_glc_e']),
        'biomass': abs(res_fba['UOK262']['biomass'].pre_solution.fobj / res_env['UOK262']['R_EX_glc_e'])
    }
])

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.factorplot('condition', 'atp', data=plot_df, palette=pal, kind='bar', lw=0)
plt.axhline(4, ls='--', lw=.3, c='gray')
plt.axhline(32, ls='--', lw=.3, c='gray')
plt.gcf().set_size_inches(1, 2)
plt.ylabel('ATP yield (per mol Glucose)')
plt.xlabel('')
plt.savefig('./reports/atp_yield_bar.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.factorplot('condition', 'biomass', data=plot_df, palette=pal, kind='bar', lw=0)
plt.gcf().set_size_inches(1, 2)
plt.ylabel('Biomass yield (per mol Glucose)')
plt.xlabel('')
plt.savefig('./reports/biomass_yield_bar.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# O2/CO2 rates
plot_df = DataFrame([{'condition': c, 'o2': res_fba[c]['atp'].values['R_EX_o2_e'], 'co2': res_fba[c]['atp'].values['R_EX_co2_e']} for c in conditions])
plot_df['o2'] *= -1

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.regplot('co2', 'o2', plot_df[plot_df['condition'] == 'UOK262'], fit_reg=False, color=pal['UOK262'], label='UOK262')
sns.regplot('co2', 'o2', plot_df[plot_df['condition'] == 'UOK262pFH'], fit_reg=False, color=pal['UOK262pFH'], label='UOK262pFH')
sns.despine()
plt.xlabel('CO2 secretion (mmol/gDW/h)')
plt.ylabel('O2 consumption (mmol/gDW/h)')
plt.legend(loc=3)
plt.gcf().set_size_inches(2, 2)
plt.savefig('./reports/ocr_ecar_scatter.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Export simulations
fluxes = DataFrame({c: res_fba[c]['atp'].values for c in res_fba})
fluxes.to_csv('./data/pfba_atp.csv')
