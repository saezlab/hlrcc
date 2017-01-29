#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from pandas import Series, DataFrame

# -- Import label-free MS
proteomics_lf = Series.from_csv('./data/uok262_proteomics_labelfree_processed.csv')
phosphoproteomics_lf = Series.from_csv('./data/uok262_phosphoproteomics_labelfree_processed.csv')

# -- Import TMT-MS
proteomics_tmt = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')
phosphoproteomics_tmt = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')


# -- Correlate: phosphoproteomics
plot_df = DataFrame({'label-free': phosphoproteomics_lf, 'tmt': phosphoproteomics_tmt}).dropna()

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
g.set_axis_labels('Label-free (log2 FC)', 'TMT (log2 FC)')
plt.savefig('./reports/phosphoproteomics_tmt_labelfree_correlation.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'


# -- Correlate: proteomics
plot_df = DataFrame({'label-free': proteomics_lf, 'tmt': proteomics_tmt}).dropna()

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
g.set_axis_labels('Label-free (log2 FC)', 'TMT (log2 FC)')
plt.savefig('./reports/proteomics_tmt_labelfree_correlation.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
