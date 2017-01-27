#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import rpathways
from scipy.stats.stats import pearsonr
from pandas import read_csv, DataFrame, Series
from hlrcc.metabolomics.sampler import sample
from framed import load_cbmodel, simplify, FBA, MOMA, pFBA, FVA, lMOMA, reaction_deletion, gene_deletion, essential_genes
from framed.cobra.deletion import deleted_genes_to_reactions

# Rna-seq
genename_map = read_csv('./files/protein-coding_gene.txt', sep='\t')
genename_map['hgsn'] = ['G_' + i.replace(':', '_') for i in genename_map['hgnc_id']]
genename_map = genename_map.groupby('hgsn')['symbol'].agg(lambda x: list(x))

gmap = read_csv('./files/protein-coding_gene.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('ensembl_gene_id')['hgsn'].agg(lambda x: list(x))

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

# Not-expressed transcripts
transcripts = read_csv('./data/UOK262_expressed_transcripts.csv')['genes']
transcripts = [gmap[i][0] for i in transcripts if i in gmap]
inactive_genes = {g for g in model.genes if g not in transcripts}
inactive_reactions = deleted_genes_to_reactions(model, inactive_genes)

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
print '[INFO] Exported'


# --
fluxes = DataFrame({c: res_fba[c]['atp'].values for c in res_fba})
fluxes['delta'] = fluxes['UOK262'].abs() - fluxes['UOK262pFH'].abs()
# fluxes = fluxes[fluxes['delta'].abs() > 1e-5]
print fluxes.sort('delta')


# -- Flux pathways heatmaps
plot_df = fluxes.ix[rpathways['glycolysis'] + rpathways['mitochondria']].dropna()

pal, cmap = dict(zip(*(['glycolysis', 'mitochondria'], sns.color_palette('Set2', n_colors=6).as_hex()))), sns.diverging_palette(10, 220, n=7, as_cmap=True)
rcol = Series({i: pal[p] for i in plot_df.index for p in rpathways if i in rpathways[p]}, name='Pathway')

sns.set(style='white', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False, figsize=(1, 6), row_colors=rcol, square=True)
plt.savefig('./reports/flux_heatmap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# # --
# gmap = read_csv('./files/protein-coding_gene.txt', sep='\t')
# gmap['hgnc_id'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
# gmap = gmap.set_index('hgnc_id')['symbol']
#
# umap = read_csv('./files/protein-coding_gene.txt', sep='\t').dropna(subset=['uniprot_ids'])
# umap = umap.groupby('uniprot_ids')['symbol'].agg(lambda x: ';'.join([g for i in x for g in i.split('|')]))
#
# # # Proteomics
# # proteomics = read_csv('./data/uok262_proteomics_logfc.txt', sep='\t')
# # proteomics['gene'] = [umap.ix[i] if i in umap.index else np.nan for i in proteomics.index]
# # proteomics_fc = proteomics.groupby('gene')['logFC'].mean()
# # print proteomics_fc.sort_values()
#
# # Proteomics
# proteomics = read_csv('./data/uok262_proteomics_tmt.csv')
# proteomics['gene'] = [umap.ix[i] if i in umap.index else np.nan for i in proteomics['Accession']]
# proteomics = proteomics.dropna()
# proteomics = proteomics.groupby('gene').median().mean(1)
# proteomics_fc = proteomics.copy()
# print proteomics_fc.sort_values()
#
# # # Phosphoproteomics
# # phospho = read_csv('./data/uok262_phosphoproteomics_logfc.txt', sep='\t')
# # phospho['psite'] = ['%s_%s' % (umap[i.split('_')[0]], i.split('_')[1]) if i.split('_')[0] in umap else np.nan for i in phospho.index]
# # phospho = phospho.dropna()
# # phospho['gene'] = [i.split('_')[0] for i in phospho['psite']]
# # phospho_fc = phospho.groupby('psite')['logFC'].median()
# # print phospho_fc.sort_values()
# #
# # psites_map = phospho.groupby('gene')['psite'].agg(lambda x: set(x)).to_dict()
#
# r_genes = {r: {gmap.ix[g] for g in model.reactions[r].gpr.get_genes() if g in gmap.index} for r in model.reactions if model.reactions[r].gpr}
# pathways = {r: model.reactions[r].metadata['SUBSYSTEM'] for r in model.reactions if 'SUBSYSTEM' in model.reactions[r].metadata and model.reactions[r].metadata['SUBSYSTEM'] != ''}
#
# plot_df = DataFrame([
#     {'reaction': r, 'gene': g, 'pathway': pathways[r], 'flux': fluxes.ix[r, 'delta'], 'protein': proteomics_fc[g]}
#     for r in fluxes.index if r in fluxes.index and r in r_genes and r in pathways for g in r_genes[r] if g in proteomics_fc
# ])
#
# plot_df = plot_df[plot_df['flux'].abs() > 1e-5]
# plot_df = plot_df.groupby('pathway').mean().sort('protein')
# print plot_df
#
# cor, pval = pearsonr(plot_df['flux'], plot_df['protein'])
# print cor, pval
#
# p_highlight = ['Oxidative phosphorylation', 'NAD metabolism', 'Citric acid cycle', 'Glutamate metabolism', 'Cysteine Metabolism', 'Glycine, serine, alanine and threonine metabolism']
# pal = dict(zip(*(p_highlight + ['Others'], sns.color_palette('Set2', n_colors=6).as_hex() + ['#D3D3D3'])))
#
# sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# sns.regplot('flux', 'protein', plot_df, fit_reg=True, color='#34495e', line_kws={'lw': .3})
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
