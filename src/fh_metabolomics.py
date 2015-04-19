__author__ = 'emanuel'

import operator
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import DataFrame, read_csv, melt
from utils.plot_utils import save_plot

# General variables
data_dir = '/Users/emanuel/Projects/data/fh_cells/'
conditions = ['UOK262', 'UOK262pFH']

# Import medium
medium = read_csv(data_dir + 'human_metabolomics/human_metabolomics.txt', sep='\t')
medium = medium.drop('exchange', axis=1)
medium.columns = [c.split('_')[0] for c in medium.columns]

# Parse table to data-frame
medium_melted = melt(medium, id_vars=['metabolite'])

# Col order
condition_dif = {met: np.abs(np.median(medium_melted.loc[np.bitwise_and(medium_melted['metabolite'] == met, medium_melted['variable'] == 'UOK262'), 'value'].values) - np.median(medium_melted.loc[np.bitwise_and(medium_melted['metabolite'] == met, medium_melted['variable'] == 'UOK262pFH'), 'value'].values)) for met in set(medium_melted['metabolite'].values)}
condition_dif = [k for k, v in sorted(condition_dif.items(), key=operator.itemgetter(1), reverse=True)]

# Plot
sns.set(style='white')
sns.set_context('paper')
g = sns.factorplot('metabolite', 'value', 'variable', estimator=np.median, data=medium_melted, kind='point', palette=sns.color_palette('Paired'), x_order=condition_dif, aspect=2, legend=False, ci=None)
g.set_ylabels('mol / gDW / h')
g.set_xticklabels(rotation=30)
plt.legend(loc='lower right')
plt.axhline(0, c='#95a5a6', lw=.5, alpha=.15)
g.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/' + 'metabolic_consumption_release' + '.pdf')