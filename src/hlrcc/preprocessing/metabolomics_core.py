import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from sklearn.linear_model import LinearRegression
from pandas import DataFrame, Series, read_csv

# Import metabolite map
m_map = read_csv('%s/files/metabolites_map.txt' % wd, sep='\t', index_col=1)
m_map.index = [i.lower() for i in m_map.index]
m_map = m_map.to_dict()['exchange']

# Import data-sets
core_1 = read_csv('%s/data/metabolomics_core_replicate_1.txt' % wd, sep='\t', index_col=0).T
core_2 = read_csv('%s/data/metabolomics_core_replicate_2.txt' % wd, sep='\t', index_col=0).T

# Overlap data-sets
core_1.index = [i.lower() for i in core_1.index]
core_2.index = [i.lower() for i in core_2.index]

core_1.columns = ['%s_%d_rep1' % (core_1.columns[i - 1], i) for i in range(1, len(core_1.columns) + 1)]
core_2.columns = ['%s_%d_rep2' % (core_2.columns[i - 1], i) for i in range(1, len(core_2.columns) + 1)]

metabolites = set(core_1.index).intersection(core_2.index)
print '[INFO] Metabolites measured: ', len(metabolites)

# Merge two data-sets
core = core_1.join(core_2) * 1000

# # Filter by variation
# core = core[core.std(1) > .1]

# Map metabolite to exchange reaction
core.index = [m_map[i] for i in core.index]

# Export data-set
core.to_csv('%s/data/uok262_metabolomics_core_processed.txt' % wd, sep='\t')
print '[INFO] Export metabolomics'

# Remove batch experiment effect
design = DataFrame({'replicate%d' % i: {c: int(c.endswith('_rep%d' % i)) for c in core} for i in range(1, 3)}).ix[core.columns]


def rm_batch(x, y):
    ys = y.dropna()
    xs = x.ix[ys.index]

    lm = LinearRegression().fit(xs, ys)

    return ys - xs.dot(lm.coef_) - lm.intercept_

core = DataFrame({metabolite: rm_batch(design, core.ix[metabolite, core.columns]) for metabolite in core.index}).T
print '[INFO] Replicate effect removed: ', core.shape

# -- Plot
# Heatmap
rep_cmap = dict(zip(*(['rep1', 'rep2'], sns.light_palette('#e74c3c', 3)[1:])))
cod_cmap = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#3498db', 3)[1:])))

plot_df = core.corr(method='spearman')

row_color = [[rep_cmap[i.split('_')[2]] for i in plot_df.index], [cod_cmap[i.split('_')[0]] for i in plot_df.index]]
col_color = [[rep_cmap[i.split('_')[2]] for i in plot_df], [cod_cmap[i.split('_')[0]] for i in plot_df]]

cmap = sns.light_palette('#34495e', 10, as_cmap=True)
sns.set(style='white', context='paper', font_scale=0.75)
sns.clustermap(plot_df, annot=True, cmap=cmap, lw=.3, row_colors=row_color, col_colors=col_color)
plt.savefig('%s/reports/metabolomics_core_clutermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'