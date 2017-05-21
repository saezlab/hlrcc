#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pandas import DataFrame, Series, read_csv


# -- Imports
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')
proteomics = read_csv('./data/uok262_proteomics_labelfree_processed_fc.csv', index_col=0)


# -- Correlate
plot_df = DataFrame({'transcriptomics': transcriptomics, 'proteomics': proteomics['fc']}).dropna()

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
g = sns.jointplot(
    'proteomics', 'transcriptomics', plot_df, 'reg', color='#34495e', space=0,
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}, 'fit_reg': True},
    marginal_kws={'hist': False, 'rug': False, 'kde': False}, stat_func=spearmanr,
    annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4}
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Proteomics (log2 FC)', 'Transcriptomics (log2 FC)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/proteomics_transcriptomics_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
