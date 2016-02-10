import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from hlrcc import wd
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv
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

    lm = sm.OLS(ys, xs).fit_regularized(L1_wt=0, alpha=alpha)

    return lm.params.drop('Const')


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

plot_df = DataFrame({'lm': k_activity_lm.ix[kinases_ov], 'gsea': k_activity_gsea.ix[kinases_ov]})

sns.set(style='ticks')
g = sns.jointplot(
    'lm', 'gsea', plot_df, 'reg', color='#34495e', joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5}},
    marginal_kws={'hist': False, 'rug': True}, annot_kws={'template': 'Pearson: {val:.2g}, p-value: {p:.1e}'}, ylim=[-5, 5], space=0
)
plt.axhline(0, ls='--', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='--', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Kinase activities (Ridge)', 'Kinase activities (GSEA)')
plt.savefig('%s/reports/kinases_activites_jointplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'


# -- MSigDB pathway enrichment
signatures = read_gmt('%s/files/c2.cp.kegg.v5.1.symbols.gmt' % wd)
kinase_set = {human_uniprot[k][0] for k, v in k_activity_lm.to_dict().items() if k in human_uniprot}
kinase_all = {human_uniprot[k][0] for k in k_targets if k in human_uniprot}

# hypergeom.sf(x, M, n, N, loc=0)
# M: total number of objects,
# n: total number of type I objects
# N: total number of type I objects drawn without replacement
pathways_es = {sig: (hypergeom.sf(
    len(kinase_set.intersection(signatures[sig])),
    len(kinase_all),
    len(kinase_all.intersection(signatures[sig])),
    len(kinase_set)
), len(kinase_set.intersection(signatures[sig]))) for sig in signatures if len(kinase_set.intersection(signatures[sig])) > 2}
print '[INFO] GO terms enrichment done'

pathways_es = DataFrame(pathways_es, index=['pvalue', 'intersection']).T.dropna()
pathways_es = pathways_es[pathways_es['pvalue'] != 0]
pathways_es['fdr'] = multipletests(pathways_es['pvalue'], method='fdr_bh')[1]
pathways_es = pathways_es.sort('fdr', ascending=False)
pathways_es = pathways_es[pathways_es['pvalue'] < 0.05]
pathways_es['name'] = [' '.join(i.lower().split('_')[1:]) for i in pathways_es.index]
print pathways_es

sns.set(style='ticks')
plot_df = pathways_es[pathways_es['pvalue'] < .05]
colours, y_pos = sns.color_palette('Paired', 2), [x + 1.5 for x in range(len(plot_df['name']))]

plt.barh(y_pos, -np.log10(plot_df['pvalue']), lw=0, align='center', height=.5, color=colours[0], label='p-value')
plt.barh(y_pos, -np.log10(plot_df['fdr']), lw=0, align='center', height=.5, color=colours[1], label='FDR')
plt.yticks(y_pos, plot_df['name'])

plt.axvline(-np.log10(0.05), ls='--', lw=0.4, c='gray')
plt.axvline(-np.log10(0.01), ls='--', lw=0.4, c='gray')

plt.text(-np.log10(0.05) * 1.01, .5, '5%', ha='left', color='gray', fontsize=9)
plt.text(-np.log10(0.01) * 1.01, .5, '1%', ha='left', color='gray', fontsize=9)

sns.despine()
plt.xlabel('-log10')
plt.title('KEGG pathways enrichment')
plt.legend(loc=4)
plt.gcf().set_size_inches(5., 8., forward=True)
plt.savefig('%s/reports/k_activity_pathway_enrichment.pdf' % wd, bbox_inches='tight')
plt.close()
print '[INFO] Pathways enrichment plotted'

