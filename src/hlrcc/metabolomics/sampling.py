#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import rpathways
from pymist.enrichment.gsea import gsea
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from scipy.stats.stats import pearsonr
from pandas import read_csv, DataFrame, Series
from hlrcc.metabolomics.sampler import sample
from framed import load_cbmodel, simplify, FBA, pFBA, FVA, lMOMA, reaction_deletion, gene_deletion, essential_genes
from framed.cobra.deletion import deleted_genes_to_reactions

# -- Gene map
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])


# -- Import metabolic model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# FVA to simplify model
simplify(model)
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))


# -- Imports
# RNA-seq
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')

# Proteomics
proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')

# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.csv', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .1]

# Medium
medium = read_csv('./files/DMEM_41966_medium_revised.txt', sep='\t').dropna().set_index('exchange')

# O2 consumption (mmol/gDW/h)
o2_exch = 'R_EX_o2_e'
core_o2 = read_csv('./data/uok262_metabolomics_core_o2_processed.txt', sep='\t', index_col=0)[o2_exch]
core_o2 *= -1


# -- Non-expressed transcripts
inactive_genes = {g for g in model.genes if gmap[g] not in transcriptomics and gmap[g] not in proteomics}
inactive_reactions = deleted_genes_to_reactions(model, inactive_genes)
print 'Inactive %d genes, %d reactions' % (len(inactive_genes), len(inactive_reactions))


# -- Fit medium
conditions = ['UOK262', 'UOK262pFH']

res_fba, res_env = {}, {}

# c = 'UOK262'
for c in conditions:
    print c
    res_fba[c] = {}

    # MOMA with environmental conditions restricted to metabolites in the medium
    env = {r: (medium.ix[r, 'lb'], model.reactions[r].ub) if r in list(medium.index) else (0, model.reactions[r].ub) for r in model.reactions if r.startswith('R_EX_') or r.startswith('R_sink_') or r.startswith('R_DM_')}

    # Add FH KO to tumour cell lines
    if c == 'UOK262':
        env.update({'R_FUM': (0, 0), 'R_FUMm': (0, 0)})

    # Fix Biomass
    env.update({'R_biomass_reaction': 0.02721 if c == 'UOK262' else 0.02486})

    # Minimise differences of measured rates: [(r, meas[r], moma_sol.values[r]) for r in meas]
    meas = {r: core.ix[r, c] if r != o2_exch else core_o2[c] for r in list(core.index) if r in model.reactions}
    moma_sol = lMOMA(model, reference=meas, constraints=env)
    env.update({r: moma_sol.values[r] for r in meas})
    print moma_sol

    # Save environmental conditions
    res_env[c] = env

    # Biomass and ATP production
    res_fba[c]['biomass'] = pFBA(model, objective={'R_biomass_reaction': 1}, constraints=env)
    print 'Max biomass: ', res_fba[c]['biomass'].pre_solution.fobj

    res_fba[c]['atp'] = pFBA(model, objective={'R_ATPM': 1}, constraints=env)
    print 'Max ATP: ', res_fba[c]['atp'].pre_solution.fobj

    # solution = res_fba[c]['atp']
    # print solution.show_metabolite_balance('M_co2_m', model)
    # print solution.values['R_CO2tm']

    # # Sampler
    # sampling = sample(model, n_samples=200, n_steps=1000, verbose=1, constraints=env)
    # sampling.to_csv('./data/%s_sampling.csv' % c, index=False)
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
plt.yticks(range(0, 36, 4))
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


# -- Store fluxes
fluxes = DataFrame({c: res_fba[c]['atp'].values for c in res_fba})
fluxes['delta'] = fluxes['UOK262'].abs() - fluxes['UOK262pFH'].abs()
fluxes.to_csv('./data/pfba_atp.csv')


# -- Fluxes from sampling
ko_sampling, wt_sampling = [read_csv('./data/%s_sampling.txt' % c, sep='\t') for c in ['UOK262', 'UOK262pFH']]

fluxes = DataFrame(
    {r: {'UOK262': ko_sampling[r].median(), 'UOK262pFH': wt_sampling[r].median()} for r in ko_sampling}
).T
fluxes['delta'] = fluxes['UOK262'].abs() - fluxes['UOK262pFH'].abs()
print fluxes.sort('delta')


# -- Fluxes heatmaps
cmap = sns.diverging_palette(10, 220, sep=5, n=20, as_cmap=True)

# Flux pathway heatmap
plot_df = fluxes.ix[rpathways['glycolysis'] + rpathways['mitochondria']].dropna()
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

pal = dict(zip(*(['glycolysis', 'mitochondria'], sns.color_palette('Set2', n_colors=6).as_hex())))
rcol = Series({i: pal[p] for i in plot_df.index for p in rpathways if i in rpathways[p]}, name='')

sns.set(style='white', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False, row_colors=rcol)
plt.gcf().set_size_inches(.5, 5)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'

# Exchange reactions
plot_df = fluxes.ix[{r for c in conditions for r in res_env[c] if res_fba[c]['atp'].values[r] != 0}].dropna().sort_values('UOK262')
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False)
plt.gcf().set_size_inches(.5, 4)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap_exchange.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# -- Protein vs Flux
# Maps
r_genes = {r: {gmap.ix[g] for g in model.reactions[r].gpr.get_genes() if g in gmap.index} for r in model.reactions if model.reactions[r].gpr}

p_reactions = DataFrame([{'r': r, 'p': model.reactions[r].metadata['SUBSYSTEM']} for r in model.reactions if 'SUBSYSTEM' in model.reactions[r].metadata and model.reactions[r].metadata['SUBSYSTEM'] != ''])
p_reactions = p_reactions.groupby('p')['r'].agg(lambda x: set(x)).to_dict()

p_reactions_genes = {p: {g for r in p_reactions[p] if r in r_genes for g in r_genes[r]} for p in p_reactions}

# Enrichment pathway analysis - Fluxomics
dataset = fluxes['delta'].replace(0, np.nan).dropna().to_dict()

df_enrichment = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, 1)) for c, sig in p_reactions.items()]
df_enrichment = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment]).dropna().set_index('pathway')
df_enrichment.sort(['escore']).to_csv('./files/fluxomics_pathways.csv', index=False)
# df_enrichment = read_csv('./files/fluxomics_pathways.csv')

# Enrichment pathway analysis - Proteomics
dataset = proteomics.to_dict()

df_enrichment_p = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, 1000)) for c, sig in p_reactions_genes.items()]
df_enrichment_p = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment_p]).dropna().set_index('pathway')
df_enrichment_p.sort(['escore']).to_csv('./files/proteomics_pathways.csv', index=False)
# df_enrichment = read_csv('./files/proteomics_pathways.csv')

p_overlap = set(df_enrichment[df_enrichment['length'] > 2].index).intersection(df_enrichment_p[df_enrichment_p['length'] > 2].index)
print df_enrichment.ix[p_overlap].sort_values(['escore'])
print df_enrichment_p.ix[p_overlap].sort_values(['escore'])

print df_enrichment.ix[['Heme synthesis', 'Heme degradation']]


# #
# plot_df = DataFrame([
#     {
#         'flux': df_enrichment.ix[p, 'escore'],
#         'protein': df_enrichment_p.ix[p, 'escore'],
#         'pathway': p
#     } for p in set(df_enrichment_p.index).intersection(df_enrichment.index)
# ]).dropna().set_index('pathway')
# plot_df = plot_df[plot_df['flux'].abs() > 1e-5]
# print plot_df.sort('flux')
#
# cor, pval = pearsonr(plot_df['flux'], plot_df['protein'])
# print cor, pval
#
# p_highlight = [
#     'Oxidative phosphorylation',
#     'NAD metabolism',
#     'Citric acid cycle',
#     'Glutamate metabolism',
#     'Cysteine Metabolism',
#     'Glycine, serine, alanine and threonine metabolism'
# ]
# pal = dict(zip(*(p_highlight + ['Others'], sns.color_palette('Set2', n_colors=6).as_hex() + ['#D3D3D3'])))
#
# # Plot
# sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# sns.regplot('flux', 'protein', plot_df, fit_reg=False, color='#34495e', line_kws={'lw': .3})
# for p in plot_df.index:
#     if p in p_highlight:
#         sns.regplot('flux', 'protein', plot_df.ix[[p]], fit_reg=False, label=p, color=pal[p], scatter_kws={'s': 40, 'lw': 0})
#     else:
#         sns.regplot('flux', 'protein', plot_df.ix[[p]], fit_reg=False, color=pal['Others'], scatter_kws={'s': 40, 'lw': 0})
#
# sns.despine()
# plt.axhline(0, ls='-', lw=.3, c='gray')
# plt.axvline(0, ls='-', lw=.3, c='gray')
# plt.xlabel('Flux mmol/gDW/h (UOK262 - UOK262pFH)')
# plt.ylabel('Protein log2 FC (UOK262 - UOK262pFH)')
# plt.legend(loc=3)
# plt.title('Pearson r: %.2f, p-value: %.2e' % (cor, pval))
# plt.gcf().set_size_inches(3, 3)
# plt.savefig('./reports/protein_flux_scatter.pdf', bbox_inches='tight')
# plt.close('all')
# print '[INFO] Plot done'


# sns.set(style='ticks', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.direction': 'out', 'ytick.direction': 'out'})
# g = sns.factorplot(x='escore', y='pathway', data=df_enrichment, kind='bar', ci=None, legend_out=False, lw=.3, orient='h', color='#99A3A4')
# g.despine(trim=True)
# plt.xlabel('Enrichment score')
# plt.ylabel('')
# plt.gcf().set_size_inches(1, 2)
# # plt.xticks(np.arange(.0, .3, .1))
# plt.title('Metabolic pathways (KO vs WT)')
# plt.savefig('./reports/fluxomics_pathways_enrichment_barplot.pdf', bbox_inches='tight')
# plt.close('all')
# print '[INFO] Fluxomics'