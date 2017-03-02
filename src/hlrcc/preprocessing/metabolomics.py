#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from scipy.stats.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model.base import LinearRegression

# ---- Quantify metabolite consumption/release (umol/ugDW/h)
# -- Import protein dry weights (ugDW/well/h) - 75% DW is protein
dry_weights = read_csv('./data/core/Data_CoRe4_protein_dw.csv', index_col=0)

# -- Build linear regrssions for metabolites standards concentration
medium_s = read_csv('./data/core/Data_CoRe4_all_standards.csv', index_col=0)
metabo_s = read_csv('./data/core/Data_CoRe4_all_standards_metabolites.csv').dropna()

metabolites = list(set(metabo_s['metabolite']))

lm_metabolites = {}
# m = 'phenylalanine'
for m in metabolites:
    # Metabolite standard log10 a.u.
    x = metabo_s[metabo_s['metabolite'] == m].set_index('standard').drop('metabolite', axis=1)[['log10']]

    # Metabolite standard log10 concentration
    y = medium_s.ix[x.index, 'log10']

    # Fit linear regression
    lm = LinearRegression().fit(x, y)
    lm_metabolites[m] = lm
    print '%s: y = %.2fx + %.2f' % (m, lm.coef_[0], lm.intercept_)

# -- Quantify CORE raw spectra
# Import CORE raw spectra
core = read_csv('./data/core/Data_CoRe4_all.csv', index_col=1).replace(0, np.nan).drop(['HEPES'], axis=1)

# Exclude replicate 3
core = core[core['replicate'] != 3]

# Log10 transform
core[metabolites] = np.log10(core[metabolites]).replace(np.nan, 0)

# Fit concentration regression: a.u. -> uM. Note: NaNs are kept
nans = core.isnull()

for m in metabolites:
    core[[m]] = lm_metabolites[m].predict(core[[m]].replace(np.nan, 0))

core[nans] = np.nan

# Concentration in the extract: Revert log10 (uM)
core[metabolites] = 10 ** core[metabolites]

# Concentration in the well: 16 dilution (uM)
core[metabolites] *= 16

# Amount in well (u moles)
core[metabolites] = core[metabolites] * 1.5 / 1000

# Delta to control (u moles)
core = DataFrame({m: {i: core.ix[i, m] - core.ix[[x for x in core.index if x.startswith(i.split('_')[0] + '_M_')], m].mean() for i in core.index if '_M_' not in i} for m in metabolites})

# Quantification by DW (divide by ug dry weight generated/well): umol/ugDW
core = DataFrame({m: {i: core.ix[i, m] / dry_weights.ix[int(i.split('_')[0]), i.split('_')[1]] for i in core.index if '_M_' not in i} for m in metabolites})

# Divide per hour: umol/ugDW/h
core /= 24

# Convert to m mol: mmol/gDW/h
core *= 1000

# Export
core.to_csv('./data/core/Data_CoRe4_all_mmol_gDW_h.csv')
# core = read_csv('./data/core/Data_CoRe4_all_mmol_gDW_h.csv', index_col=0)
print core

# ---- Statistical analysis
core = core.reset_index()
core['condition'] = ['UOK262pFH' if i.split('_')[1] == 'pFH' else 'UOK262' for i in core['index']]
core = core.set_index(['condition', 'index'])

# -- Metabolite exchange reaction map
m_map = read_csv('./files/metabolites_map.txt', sep='\t', index_col=1)
m_map = m_map.to_dict()['exchange']

# -- Differential rates between conditions
core_ttest = {m: ttest_ind(core.ix['UOK262'][m], core.ix['UOK262pFH'][m], equal_var=False) for m in metabolites}
core_ttest = DataFrame([{
    'metabolite': m,
    't': core_ttest[m][0],
    'pval': core_ttest[m][1],
    'diff': core.ix['UOK262'][m].median() - core.ix['UOK262pFH'][m].median(),
    'UOK262': core.ix['UOK262'][m].median(),
    'UOK262pFH': core.ix['UOK262pFH'][m].median(),
} for m in core]).dropna()
core_ttest['fdr'] = multipletests(core_ttest['pval'], method='fdr_bh')[1]
core_ttest['exchange'] = [m_map[m] if m in m_map else np.nan for m in core_ttest['metabolite']]
print core_ttest.sort('fdr')

# -- Boxplots
fdr_thres = 0.05

plot_df = core[core_ttest[core_ttest['fdr'] < fdr_thres]['metabolite']].reset_index().drop('index', axis=1).set_index('condition').unstack().reset_index()
plot_df.columns = ['metabolite', 'condition', 'rate']
plot_df = plot_df.sort_values('rate', ascending=False)
plot_df['metabolite'] = [i.capitalize() for i in plot_df['metabolite']]

order = [i.capitalize() for i in core_ttest[core_ttest['fdr'] < fdr_thres].sort('diff', ascending=False)['metabolite']]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, col='metabolite', col_order=order, col_wrap=7, size=1.85, aspect=.5, legend_out=True, sharex=False, sharey=False)
g.map(sns.boxplot, 'condition', 'rate', orient='v', palette=sns.light_palette('#34495e', 3)[1:], linewidth=.3, order=['UOK262', 'UOK262pFH'], fliersize=2)
g.map(plt.axhline, y=0, ls='-', lw=.3, alpha=.7, color='gray')
g.despine(trim=True, bottom=True)
g.set_titles('{col_name}')
g.set_xlabels('')
g.set_ylabels('Flux rate (mmol/gDW/h)')
g.set(xticks=[])
g.fig.subplots_adjust(wspace=1.5)
plt.savefig('./reports/metabolomics_core_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# -- Corr heatmap
plot_df = core[core_ttest[core_ttest['fdr'] < fdr_thres]['metabolite']].reset_index().drop('condition', axis=1).set_index('index').copy()
plot_df = plot_df.T.corr(method='spearman')

rep_cmap = dict(zip(*(['180716', '220716', '280716'], sns.light_palette('#e74c3c', 4).as_hex()[1:])))
cod_cmap = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#3498db', 3).as_hex()[1:])))
cod_map = {'262': 'UOK262', 'pFH': 'UOK262pFH'}

row_color = DataFrame({'Replicate': {i: rep_cmap[i.split('_')[0]] for i in plot_df.index}, 'Condition': {i: cod_cmap[cod_map[i.split('_')[1]]] for i in plot_df.index}})
col_color = DataFrame({'Replicate': {i: rep_cmap[i.split('_')[0]] for i in plot_df}, 'Condition': {i: cod_cmap[cod_map[i.split('_')[1]]] for i in plot_df}})

cmap = sns.light_palette('#34495e', 10, as_cmap=True)
sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(plot_df, cmap=cmap, lw=.3, row_colors=row_color, col_colors=col_color, figsize=(5, 5))
plt.savefig('./reports/metabolomics_core_clutermap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'

# -- Heatmap
plot_df = core[core_ttest[core_ttest['fdr'] < fdr_thres]['metabolite']].reset_index().drop('condition', axis=1).set_index('index').copy().T

rep_cmap = dict(zip(*(['180716', '220716', '280716'], sns.light_palette('#e74c3c', 4).as_hex()[1:])))
cod_cmap = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#3498db', 3).as_hex()[1:])))
cod_map = {'262': 'UOK262', 'pFH': 'UOK262pFH'}

col_color = DataFrame({'Replicate': {i: rep_cmap[i.split('_')[0]] for i in plot_df}, 'Condition': {i: cod_cmap[cod_map[i.split('_')[1]]] for i in plot_df}})

cmap = sns.light_palette('#34495e', 10, as_cmap=True)
sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(plot_df, cmap=cmap, lw=.3, col_colors=col_color, figsize=(5, 5), z_score=0)
plt.savefig('./reports/metabolomics_core_heatmap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'

# -- Export data-set
core_ttest.to_csv('./data/uok262_metabolomics_core_processed.csv')
print '[INFO] Export metabolomics'
