#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, Series
from scipy.stats.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# -- Metabolite exchange reaction map
m_map = read_csv('./files/metabolites_map.txt', sep='\t', index_col=1)
m_map = m_map.to_dict()['exchange']

# -- Import
metabolomics = read_csv('./data/core/Data_CoRe4.csv')
metabolites = list(metabolomics.drop(['replicate', 'sample'], axis=1))

# -- Differential rates between conditions
core = {m: ttest_ind(metabolomics.loc[metabolomics['sample'] == 'UOK262', m], metabolomics.loc[metabolomics['sample'] == 'UOK262pFH', m], equal_var=False) for m in metabolites}
core = DataFrame([{
    'metabolite': m,
    't': core[m][0],
    'pval': core[m][1],
    'median_diff': metabolomics.loc[metabolomics['sample'] == 'UOK262', m].median() - metabolomics.loc[metabolomics['sample'] == 'UOK262pFH', m].median(),
    'UOK262': metabolomics.loc[metabolomics['sample'] == 'UOK262', m].median(),
    'UOK262pFH': metabolomics.loc[metabolomics['sample'] == 'UOK262pFH', m].median(),
} for m in core])
core['fdr'] = multipletests(core['pval'], method='fdr_bh')[1]
core['exchange'] = [m_map[m] if m in m_map else np.nan for m in core['metabolite']]
print core.sort('fdr')

# --
plot_df = metabolomics.set_index('sample')[core[core['fdr'] < .05]['metabolite']].unstack().reset_index()
plot_df.columns = ['metabolite', 'condition', 'rate']

order = list(core[core['fdr'] < .05].sort('median_diff', ascending=False)['metabolite'])

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, row='metabolite', row_order=order, size=1., aspect=1., legend_out=True, sharex=False)
g.map(sns.violinplot, 'rate', 'condition', palette=sns.light_palette('#34495e', 3)[1:], split=False)
g.map(plt.axvline, x=0, ls='-', lw=.5, alpha=.7, color='gray')
g.despine(trim=True)
g.set_titles('{row_name}')
g.set_ylabels('')
g.set_xlabels('Flux rate (umol/ugDW/h)')
plt.savefig('./reports/metabolomics_core_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

# --
plot_df = metabolomics.copy()
plot_df = plot_df.set_index(plot_df['sample'] + '_' + plot_df['replicate'].astype(str) + '_' + plot_df.index.astype(str)).drop(['sample', 'replicate'], axis=1)[order].T
plot_df = plot_df.corr(method='spearman')

rep_cmap = dict(zip(*(['1', '2', '3'], sns.light_palette('#e74c3c', 4).as_hex()[1:])))
cod_cmap = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#3498db', 3).as_hex()[1:])))

row_color = DataFrame({'Replicate': {i: rep_cmap[i.split('_')[1]] for i in plot_df.index}, 'Condition': {i: cod_cmap[i.split('_')[0]] for i in plot_df.index}})
col_color = DataFrame({'Replicate': {i: rep_cmap[i.split('_')[1]] for i in plot_df}, 'Condition': {i: cod_cmap[i.split('_')[0]] for i in plot_df}})

cmap = sns.light_palette('#34495e', 10, as_cmap=True)
sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.clustermap(plot_df, cmap=cmap, lw=.3, row_colors=row_color, col_colors=col_color, figsize=(5, 5))
plt.savefig('./reports/metabolomics_core_clutermap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'


# -- Export data-set
core.to_csv('./data/uok262_metabolomics_core_processed.txt', sep='\t')
print '[INFO] Export metabolomics'

