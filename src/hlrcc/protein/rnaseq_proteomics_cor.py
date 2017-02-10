#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series


# -- Imports
proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')

# -- Correlate
ov_meas = set(transcriptomics.index).intersection(proteomics.index)

plot_df = DataFrame({'transcriptomics': transcriptomics[ov_meas], 'proteomics': proteomics[ov_meas]})

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'transcriptomics', 'proteomics', plot_df, 'reg', color='#34495e', joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}},
    marginal_kws={'hist': False, 'rug': False}, annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4}, space=0
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Transcriptomics (log2 FC)', 'Proteomics (log2 FC)')
plt.savefig('./reports/rnaseq_proteomics_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
