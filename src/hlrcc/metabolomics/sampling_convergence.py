#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from pandas import read_csv, DataFrame, Series
from statsmodels.stats.multitest import multipletests

# -- Imports
conditions = ['UOK262', 'UOK262pFH']

# Sample 1
ko_sampling_1, wt_sampling_1 = [read_csv('./data/%s_sampling_1.csv' % c) for c in ['UOK262', 'UOK262pFH']]
sampling_1 = DataFrame(
    {r: {
        'UOK262': ko_sampling_1[r].mean(),
        'UOK262pFH': wt_sampling_1[r].mean()
    } for r in set(ko_sampling_1).intersection(wt_sampling_1)}
).T
sampling_1['delta'] = sampling_1['UOK262'] - sampling_1['UOK262pFH']
print sampling_1.sort(['delta'])

# Sample 2
ko_sampling_2, wt_sampling_2 = [read_csv('./data/%s_sampling_2.csv' % c) for c in ['UOK262', 'UOK262pFH']]
sampling_2 = DataFrame(
    {r: {
        'UOK262': ko_sampling_2[r].mean(),
        'UOK262pFH': wt_sampling_2[r].mean()
    } for r in set(ko_sampling_2).intersection(wt_sampling_2)}
).T
sampling_2['delta'] = sampling_2['UOK262'] - sampling_2['UOK262pFH']
print sampling_2.sort(['delta'])


# --
reactions = set(ko_sampling_1).intersection(wt_sampling_1).intersection(ko_sampling_2).intersection(wt_sampling_2)

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(sampling_1.ix[reactions]['delta'], sampling_2.ix[reactions]['delta'])
sns.despine(trim=True)
plt.savefig('./reports/sampling_correlation.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Export list of convergent reactions
diff = (sampling_1.ix[reactions]['delta'] - sampling_2.ix[reactions]['delta']).sort_values()
diff.to_csv('./data/reactions_convergence.csv')
print diff


# sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
# plt.hist(ko_sampling_1['R_MMSAD3m'], color='r', label='UOK262', bins=30, alpha=.2)
# plt.hist(wt_sampling_1['R_MMSAD3m'], color='g', label='UOK262pFH', bins=30, alpha=.2)
# sns.despine()
# plt.legend()
# plt.xlabel('flux (mmol/gDW/h)')
# plt.title('ALDH6A1 (reaction MMSAD3m)')
# plt.gcf().set_size_inches(3, 2)
# plt.savefig('./reports/sampling_distributions.pdf', bbox_inches='tight')
# plt.close('all')
# print '[INFO] Plot done'
