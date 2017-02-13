#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, Series
from framed.solvers import solver_instance
from framed import load_cbmodel, simplify, pFBA, lMOMA, gene_deletion


# -- Imports
# Gene map
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])

# Metabolic model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Omics
proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')

phosphoproteomics = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')

transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')


# --
genes = DataFrame({g: {
    'gene': gmap[g],
    'proteomics': proteomics[gmap[g]] if gmap[g] in proteomics else np.nan,
    'transcriptomics': transcriptomics[gmap[g]] if gmap[g] in transcriptomics else np.nan
} for g in model.genes}).T
genes = genes[genes.count(1) > 1]
print genes.sort_values('proteomics')


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
res_fba, ess = {}, {}

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

    # Biomass and ATP production
    res_fba[c]['biomass'] = pFBA(c_model, objective={'R_biomass_reaction': 1})
    print 'Max biomass: ', res_fba[c]['biomass'].pre_solution.fobj

    res_fba[c]['atp'] = pFBA(c_model, objective={'R_ATPM': 1})
    print 'Max ATP: ', res_fba[c]['atp'].pre_solution.fobj

    # #
    # solution = res_fba[c]['atp']
    # print solution.show_metabolite_balance('M_3hpp_m', c_model)
    # print solution.values['R_MMMm']

    # Gene essentiallity
    # genes_2_ko = set(genes[(genes[['proteomics', 'transcriptomics']].abs() > 1.).any(1)].index)
    genes_2_ko = set(genes.index)

    c_model.set_objective({'R_ATPM': 1, 'R_biomass_reaction': 0})
    solver = solver_instance(c_model)
    ess[c] = {g: gene_deletion(c_model, [g], solver=solver) for g in genes_2_ko}
    
# ATP
plot_df = DataFrame({c: {g: ess[c][g].fobj for g in ess[c] if ess[c][g] is not None} for c in ess})
plot_df['gene'] = [gmap[g] for g in plot_df.index]
plot_df.to_csv('./data/ess_scores.csv')
print plot_df[(plot_df['UOK262pFH'] - res_fba['UOK262pFH']['atp'].pre_solution.fobj).abs() < 1e-4].sort_values('UOK262')
# print plot_df[(plot_df['UOK262'] - res_fba['UOK262']['atp'].pre_solution.fobj).abs() < 1e-4].sort_values('UOK262pFH')


# -- Plot
# Essentiallity scatter
p_highlight = ['G_HGNC_4908', 'G_HGNC_7179', 'G_HGNC_16732', 'G_HGNC_94', 'G_HGNC_3481', 'G_HGNC_3482', 'G_HGNC_3700']
pal = dict(zip(*(p_highlight + ['Others'], sns.color_palette('Set2', n_colors=7).as_hex() + ['#dfdfdf'])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
plt.axhline(res_fba['UOK262pFH']['atp'].pre_solution.fobj, ls='-', lw=.9, c='#95a5a6', alpha=.5)
plt.axvline(res_fba['UOK262']['atp'].pre_solution.fobj, ls='-', lw=.9, c='#95a5a6', alpha=.5)

for p in p_highlight:
    sns.regplot('UOK262', 'UOK262pFH', plot_df.ix[[p]], fit_reg=False, label=plot_df.ix[p, 'gene'], color=pal[p], scatter_kws={'s': 15, 'lw': 1.})

for p in plot_df.index:
    if p not in p_highlight:
        sns.regplot('UOK262', 'UOK262pFH', plot_df.ix[[p]], fit_reg=False, color=pal['Others'], scatter_kws={'s': 10, 'lw': .0, 'edgecolors': 'none'})

sns.despine(trim=True)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Enzyme')
plt.xlabel('KO (ATP production)')
plt.ylabel('WT (ATP production)')
plt.gcf().set_size_inches(2, 2)
plt.savefig('./reports/gene_lethal_scatter_atp.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# Fluxes barplot
genes = ['G_HGNC_4908', 'G_HGNC_7179', 'G_HGNC_16732']

conditions = ['KO', 'WT']
conditions_map = {'UOK262': 'KO', 'UOK262pFH': 'WT'}

fluxes = read_csv('./data/pfba_atp.csv', index_col=0).replace(np.nan, 0)
fluxes['delta'] = fluxes['UOK262'] - fluxes['UOK262pFH']
fluxes.columns = [conditions_map[c] if c in conditions_map else c for c in fluxes]

df = DataFrame([{'gene': gmap[g], 'reaction': r, 'condition': c, 'flux': fluxes.ix[r, c]} for g in genes for r in model.get_reactions_by_gene(g) if r in fluxes.index for c in conditions])
df = df[df['reaction'] != 'R_r0779']

pal = dict(zip(*(conditions, sns.light_palette('#34495e', 3, reverse=True).as_hex()[:-1])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(data=df, sharex=False)
g.map(sns.barplot, 'reaction', 'flux', 'condition', palette=pal, ci=False, lw=0)
g.set_titles('{col_name}')
g.set_xlabels('')
g.set_ylabels('Flux (mmol/gDW/h)')
# plt.legend(loc=2, title='')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('./reports/gene_lethal_flux_barplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# Sampling histograms
# Sample 1
ko_sampling_1, wt_sampling_1 = [read_csv('./data/%s_sampling_1.csv' % c) for c in ['UOK262', 'UOK262pFH']]
sampling_1 = DataFrame(
    {r: {
        'UOK262': ko_sampling_1[r].mean(),
        'UOK262pFH': wt_sampling_1[r].mean()
    } for r in set(ko_sampling_1).intersection(wt_sampling_1)}
).T
sampling_1['delta'] = sampling_1['UOK262'] - sampling_1['UOK262pFH']
print sampling_1.sort(['delta'])

# Sample 2
ko_sampling_2, wt_sampling_2 = [read_csv('./data/%s_sampling_2.csv' % c) for c in ['UOK262', 'UOK262pFH']]
sampling_2 = DataFrame(
    {r: {
        'UOK262': ko_sampling_2[r].mean(),
        'UOK262pFH': wt_sampling_2[r].mean()
    } for r in set(ko_sampling_2).intersection(wt_sampling_2)}
).T
sampling_2['delta'] = sampling_2['UOK262'] - sampling_2['UOK262pFH']
print sampling_2.sort(['delta'])

ko_sampling = ko_sampling_1.append(ko_sampling_2)
wt_sampling = wt_sampling_1.append(wt_sampling_2)

sampling = dict(zip(*(conditions, [ko_sampling, wt_sampling])))

#
df = DataFrame([{'gene': gmap[g], 'reaction': r, 'condition': c, 'flux': f} for g in genes for r in model.get_reactions_by_gene(g) for c in sampling if r in sampling[c] for f in sampling[c][r]])
df = df[df['reaction'] != 'R_r0779']

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(df, row='gene', size=.6, aspect=1, legend_out=True, sharex=True, sharey=False)
g.map(sns.boxplot, 'flux', 'reaction', 'condition', orient='h', palette=pal, linewidth=1., notch=True, fliersize=1.5)
g.map(plt.axvline, x=0, ls='-', lw=.3, alpha=.7, color='gray')
g.despine(trim=True)
g.set_ylabels('')
g.set_xlabels('Flux rate (mmol/gDW/h)')
g.set_titles('{row_name}')
plt.savefig('./reports/gene_lethal_sampling_violinplots.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
