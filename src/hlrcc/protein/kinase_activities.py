#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from pandas import Series, DataFrame
from statsmodels.stats.weightstats import ztest
from hlrcc.utils import get_ktargets, get_ktargets_omnipath


# -- Import data-sets
phospho_lf = Series.from_csv('./data/uok262_phosphoproteomics_labelfree_processed.csv')
prot_lf = Series.from_csv('./data/uok262_proteomics_labelfree_processed.csv')

phospho_tmt = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')
prot_tmt = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')


# -- Regress-out protein abundance
df_lf = DataFrame({p: {'prot': prot_lf[p.split('_')[0]], 'phospho': phospho_lf[p]} for p in phospho_lf.index if p.split('_')[0] in prot_lf}).T
lm_lf = sm.OLS(df_lf['phospho'], sm.add_constant(df_lf['prot'])).fit()
phospho_lf_res = lm_lf.resid
print lm_lf.summary()

df_tmt = DataFrame({p: {'prot': prot_tmt[p.split('_')[0]], 'phospho': phospho_tmt[p]} for p in phospho_tmt.index if p.split('_')[0] in prot_tmt}).T
lm_tmt = sm.OLS(df_tmt['phospho'], sm.add_constant(df_tmt['prot'])).fit()
phospho_tmt_res = lm_tmt.resid
print lm_tmt.summary()


# -- Import kinase-substrate interaction network
k_targets = get_ktargets_omnipath()


# -- Calculate kinase enrichment
min_targets = 1

# Kinase activity label-free
k_targets_lf = {k: k_targets[k].intersection(phospho_lf_res.index) for k in k_targets}
k_targets_lf = {k: k_targets_lf[k] for k in k_targets_lf if len(k_targets_lf[k]) >= min_targets}

k_activity_lf = DataFrame({k: dict(zip(*(['zscore', 'pvalue'], ztest(phospho_lf_res.ix[k_targets_lf[k]], phospho_lf_res.drop(k_targets_lf[k]))))) for k in k_targets_lf}).T
print k_activity_lf.sort_values('zscore')

# Kinase activity TMT-MS
k_targets_tmt = {k: k_targets[k].intersection(phospho_tmt_res.index) for k in k_targets}
k_targets_tmt = {k: k_targets_tmt[k] for k in k_targets_tmt if len(k_targets_tmt[k]) >= min_targets}

k_activity_tmt = DataFrame({k: dict(zip(*(['zscore', 'pvalue'], ztest(phospho_tmt_res.ix[k_targets_tmt[k]], phospho_tmt_res.drop(k_targets_tmt[k]))))) for k in k_targets_tmt}).T
print k_activity_tmt.sort_values('zscore')


# -- Plot kinases activities
plot_df = DataFrame({'label-free': k_activity_lf['zscore'], 'tmt': k_activity_tmt['zscore']}).dropna()

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'label-free', 'tmt', plot_df, 'reg', color='#34495e', space=0,
    marginal_kws={'hist': False, 'rug': False},
    annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4},
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}},
    xlim=[plot_df['label-free'].min() * 1.05, plot_df['label-free'].max() * 1.05],
    ylim=[plot_df['tmt'].min() * 1.05, plot_df['tmt'].max() * 1.05]
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Label-free (z-score)', 'TMT (z-score)')
plt.gcf().set_size_inches(4, 4)
plt.savefig('./reports/kinase_activities_tmt_labelfree_correlation.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'

