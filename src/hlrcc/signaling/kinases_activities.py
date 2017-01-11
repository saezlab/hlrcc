#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pymist.enrichment.gsea import gsea
from pandas.stats.misc import zscore
from pandas import DataFrame, Series, read_csv
from sklearn.linear_model.ridge import RidgeCV
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import hypergeom
from pymist.utils.read_gmt import read_gmt
from pymist.utils.omnipath_phospho import get_targets
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Kinase activity prediction function
def calc_activity_lm(x, y, alpha=0.001):
    xs = x.ix[y.index].dropna()
    xs = xs.loc[:, xs.sum() != 0]
    xs['Const'] = 1

    ys = y.ix[xs.index]

    lm = RidgeCV().fit(xs, ys)

    return zscore(Series(dict(zip(*(xs.columns, lm.coef_)))))


def calc_activity_gsea(x, y, permutations=10000):
    xs = x.replace(0, np.nan)
    xs = {i: set(xs[i].dropna().index) for i in xs}

    ys = y.to_dict()

    es = {k: gsea(ys, xs[k], permutations) for k in xs}
    es = {k: -np.log10(es[k][1]) if es[k][0] < 0 else np.log10(es[k][1]) for k in es}

    return Series(es)


# -- Import Uniprot id mapping
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)


# -- Import kinases targets
sources = ['HPRD', 'PhosphoSite', 'Signor', 'phosphoELM', 'DEPOD']
k_targets = get_targets(sources, remove_self=False)
print '[INFO] Kinases targets imported: ', k_targets.shape


# -- Import phosphoproteomics data-set
pp = read_csv('%s/data/uok262_phosphoproteomics_logfc.txt' % wd, sep='\t')
pp['psite'] = ['_'.join(i.split('_')[:2]) for i in pp.index]
pp = pp.groupby('psite')['logFC'].median()


# -- Calculate kinase enrichment
k_activity_lm = calc_activity_lm(k_targets, pp).sort(inplace=False).dropna()
k_activity_gsea = calc_activity_gsea(k_targets, pp).sort(inplace=False).dropna()

k_activity_lm.to_csv('%s/data/uok262_kinases_activity_lm.txt' % wd, sep='\t')
k_activity_gsea.to_csv('%s/data/uok262_kinases_activity_gsea.txt' % wd, sep='\t')
print '[INFO] Kinases activities estimated'

k_activity_lm = read_csv('%s/data/uok262_kinases_activity_lm.txt' % wd, sep='\t', index_col=0, names=['kinase', 'activity'])
k_activity_gsea = read_csv('%s/data/uok262_kinases_activity_gsea.txt' % wd, sep='\t', index_col=0, names=['kinase', 'activity'])

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
plt.savefig('%s/reports/kinases_activites_jointplot.pdf' % wd, bbox_inches='tight')
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
plt.savefig('%s/reports/kinases_activites_barplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
