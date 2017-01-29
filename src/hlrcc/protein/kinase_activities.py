#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.weightstats import ztest
from hlrcc.utils import get_ktargets, get_ktargets_omnipath


# -- Import data-sets
prot = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')
phospho = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')
print phospho.sort_values()

# -- Correlation with protein abundance
pp_df = DataFrame({i: {'phospho': phospho.ix[i], 'prot': prot.ix[i.split('_')[0]]} for i in phospho.index if i.split('_')[0] in prot.index}).T

lm = sm.OLS(pp_df['phospho'], sm.add_constant(pp_df['prot'])).fit()
print lm.summary()

residuals = lm.resid.copy()
print residuals.sort_values()

# -- Import kinase-substrate interaction network
k_targets = get_ktargets_omnipath()
k_targets = {k: k_targets[k].intersection(residuals.index) for k in k_targets}
k_targets = {k: k_targets[k] for k in k_targets if len(k_targets[k]) > 5}

# -- Calculate kinase enrichment
k_activity = DataFrame({k: dict(zip(*(['zscore', 'pvalue'], ztest(residuals.ix[k_targets[k]], residuals.drop(k_targets[k]))))) for k in k_targets}).T
print k_activity.sort('pvalue')


# -- Plot kinases activities
# Corrplot
kinases_ov = list(set(k_activity_lm.index).intersection(k_activity_gsea.index))

plot_df = DataFrame({'lm': k_activity_lm.ix[kinases_ov].to_dict()['activity'], 'gsea': k_activity_gsea.ix[kinases_ov].to_dict()['activity']})

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'lm', 'gsea', plot_df, 'reg', color='#34495e', joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5}},
    marginal_kws={'hist': False, 'rug': True}, annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}', 'loc': 4}, ylim=[-5, 5], space=0
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Kinase activities (Ridge)', 'Kinase activities (GSEA)')
plt.savefig('./reports/kinases_activites_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'

# Barplot
plot_df = k_activity_lm.reset_index()
plot_df = plot_df[plot_df['activity'].abs() > .5]
plot_df = plot_df[[i in human_uniprot for i in plot_df['kinase']]]
plot_df['name'] = [human_uniprot[i][0] for i in plot_df['kinase']]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.barplot('activity', 'name', data=plot_df, color='#34495e', lw=0)
plt.axvline(0, ls='-', lw=.3, alpha=.7, c='gray')
sns.despine(trim=True)
plt.xlabel('Kinase activities (GSEA)')
plt.ylabel('')
plt.title('Top kinases/phosphatases activities')
plt.gcf().set_size_inches(3., 5., forward=True)
plt.savefig('./reports/kinases_activites_barplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
