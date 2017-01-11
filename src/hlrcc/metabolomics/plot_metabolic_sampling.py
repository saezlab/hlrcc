#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pandas import DataFrame, Series, read_csv

# -- Imports
ko_sampling, wt_sampling = [read_csv('%s/data/%s_sampling.txt' % (wd, c), sep='\t', index_col=0) for c in ['UOK262', 'UOK262pFH']]


# -- Plot
reactions = ['R_HMGCOASi', 'R_GAPD', 'R_PDHm', 'R_PI35P3P', 'R_PI3PP', 'R_PGK', 'R_GAPD']

# Pairplot
ko_samples = ko_sampling[reactions]
ko_samples['condition'] = 'UOK262'

wt_samples = wt_sampling[reactions]
wt_samples['condition'] = 'UOK262pFH'

plot_df = ko_samples.append(wt_samples)

sns.set(style='ticks')
sns.pairplot(plot_df, hue='condition', palette=sns.light_palette('#34495e', 3)[1:], diag_kind='kde')
plt.savefig('%s/reports/sampling_pairplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
