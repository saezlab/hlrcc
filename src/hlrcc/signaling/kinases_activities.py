import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from hlrcc import wd
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv
from pymist.utils.omnipath_phospho import get_targets


# -- Kinase activity prediction function
def calc_activity_lm(x, y):
    xs = x.ix[y.index].dropna()
    xs['Const'] = 1

    ys = y.ix[xs.index]

    lm = sm.OLS(ys, xs).fit_regularized(L1_wt=0, alpha=0.001)

    return lm.params.drop('Const')


def calc_activity_gsea(x, y, permutations=10000):
    xs = x.replace(0, np.nan)
    xs = {i: set(xs[i].dropna().index) for i in xs}

    ys = y.to_dict()

    es = {k: gsea(ys, xs[k], permutations) for k in xs}
    es = {k: -np.log10(es[k][1]) if es[k][0] < 0 else np.log10(es[k][1]) for k in es}

    return Series(es)


# -- Import kinases targets
sources = ['HPRD', 'PhosphoSite', 'Signor', 'phosphoELM']
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


# -- Plot kinases activities
kinases_ov = list(set(k_activity_lm.index).intersection(k_activity_gsea.index))

p_highlight = ['P00519']

plot_df = DataFrame({'lm': k_activity_lm.ix[kinases_ov], 'gsea': k_activity_gsea.ix[kinases_ov]})

sns.set(style='ticks')
g = sns.jointplot(
    'lm', 'gsea', plot_df, 'reg', color='#34495e', joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5}},
    marginal_kws={'hist': False, 'rug': True}, annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}'}, ylim=[-5, 5], space=0
)
plt.axhline(0, ls='--', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='--', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
# g.ax_joint.scatter(plot_df.ix[p_highlight, 'lm'], plot_df.ix[p_highlight, 'gsea'], marker='o', c='#34495e', lw=.5, edgecolor='#e74c3c')
g.set_axis_labels('Kinase activities (Ridge)', 'Kinase activities (GSEA)')
plt.savefig('%s/reports/kinases_activites_jointplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
