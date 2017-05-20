#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from framed import load_cbmodel
from scipy.stats.stats import spearmanr
from more_itertools import unique_everseen
from pandas import read_csv, Series, DataFrame

# -- Imports
conditions = ['KO', 'WT']
conditions_map = {'UOK262': 'KO', 'UOK262pFH': 'WT'}

# Proteomics
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')
proteomics = read_csv('./data/uok262_proteomics_labelfree_processed_fc.csv', index_col=0)
phosphoproteomics = read_csv('./data/uok262_phosphoproteomics_labelfree_processed_fc.csv', index_col=0)

# Fluxomics
fluxes = read_csv('./data/pfba_atp.csv', index_col=0).replace(np.nan, 0)
fluxes['delta'] = fluxes['UOK262'] - fluxes['UOK262pFH']
fluxes.columns = [conditions_map[c] if c in conditions_map else c for c in fluxes]

# Regulatory p-sites
r_sites = read_csv('./files/Regulatory_sites.txt', sep='\t')
r_sites = r_sites[r_sites['ORGANISM'] == 'human']
r_sites = r_sites[[i.endswith('-p') for i in r_sites['MOD_RSD']]]
r_sites['res'] = ['%s_%s' % (g, p.split('-')[0]) for g, p in r_sites[['GENE', 'MOD_RSD']].values]


# -- Import metabolic model
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])

model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')

g_enzymes = {gmap[g] for g in model.genes if g in gmap}

g_to_r = model.gene_to_reaction_lookup()
g_to_r = {gmap[g]: g_to_r[g] for g in g_to_r if g in gmap}


# -- Correlation of phosphoprotoemics with proteomics
plot_df = DataFrame({p: {
    'proteomics': proteomics.loc[p.split('_')[0], 'fc'],
    'phosphoproteomics': phosphoproteomics.loc[p, 'fc']
} for p in phosphoproteomics.index if p.split('_')[0] in proteomics.index}).T

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
g = sns.jointplot(
    'proteomics', 'phosphoproteomics', plot_df, 'reg', color='#34495e', space=0,
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}, 'fit_reg': True},
    marginal_kws={'hist': False, 'rug': False, 'kde': False}, stat_func=spearmanr,
    annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4}
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Proteomics (log2 FC)', 'Phosphoproteomics (log2 FC)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/proteomics_phospho_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'

# -- Metabolic enzymes with significant phosphorylation changes for which protein abundance does not change significantly or it is not measured
plot_df = phosphoproteomics[[i.split('_')[0] in g_enzymes for i in phosphoproteomics.index]]
plot_df = plot_df[plot_df['fdr'] < 0.05]
plot_df = plot_df[[(i.split('_')[0] not in proteomics.index) or (proteomics.loc[i.split('_')[0], 'fdr'] > .05) for i in plot_df.index]]
plot_df = plot_df.reset_index().sort_values('fc')

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
sns.stripplot(y='index', x='fc', data=plot_df, color='#bfbfbf', s=4)
plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
plt.xlabel('Fold-change (log2)')
plt.ylabel('')
plt.legend(loc=0, title='Type')
plt.gcf().set_size_inches(1.5, 3.)
plt.title('Metabolic enzymes\nphosphorylation-sites')
plt.savefig('./reports/psite_regression_stripplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Match p-sites changing significantly
psites = set(plot_df['index'])

plot_df = DataFrame([{
    'Psite': p,
    'Enzyme': p.split('_')[0],
    'Reaction': r,
    'phospho (logfc)': phosphoproteomics.loc[p, 'fc'],
    'flux (delta)': fluxes.ix[r, 'delta']
 } for p in psites for r in g_to_r[p.split('_')[0]] if r in fluxes.index]).sort_values('flux (delta)')
plot_df = plot_df[plot_df['Reaction'] != 'R_ATPS4m']
plot_df = plot_df[plot_df['flux (delta)'].abs() > 1e-1]
plot_df['flux (KO)'] = [fluxes.ix[i, 'KO'] for i in plot_df['Reaction']]
plot_df['flux (WT)'] = [fluxes.ix[i, 'WT'] for i in plot_df['Reaction']]
plot_df['in_proteomics'] = [int(i in proteomics.index) for i in plot_df['Enzyme']]
plot_df['in_transcriptomics'] = [int(i in transcriptomics.index) for i in plot_df['Enzyme']]
plot_df.to_csv('./data/fluxes_phospho_table.csv', index=False)
print plot_df.sort_values('phospho (logfc)')

# # Export table
# env = Environment(loader=FileSystemLoader('.'))
# template = env.get_template('./files/table_template.html')
#
# html_out = template.render({
#     'title': 'Putative regulatory phosphorylation-sites',
#     'dataframe': plot_df.to_html()
# })
#
# HTML(string=html_out).write_pdf('./reports/regulatory_psites_table.pdf')
