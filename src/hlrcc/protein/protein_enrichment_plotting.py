#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from pandas import read_csv, DataFrame
from statsmodels.stats.multitest import multipletests


# -- Enrichment go terms
e_goterms = read_csv('./files/proteomics_tmt_go_term.csv')
e_goterms = e_goterms[e_goterms['length'] >= 3]
e_goterms['fdr'] = multipletests(e_goterms['pvalue'], method='fdr_bh')[1]
e_goterms['abs_score'] = e_goterms['escore'].abs()
print e_goterms.sort('fdr')

# Plot
plot_df = e_goterms[e_goterms['fdr'] < .05].sort_values('abs_score', ascending=False).head(20)
plot_df['name'] = [i.replace('_', ' ').lower().capitalize() for i in plot_df['signature']]
plot_df = plot_df.sort_values('escore')
print plot_df

pal = dict(zip(*(set(plot_df['type']), sns.color_palette('Set1', n_colors=3).as_hex())))

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
sns.stripplot(y='name', x='escore', hue='type', data=plot_df, palette=pal, s=4)
plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
plt.xlabel('Enrichment score\n(GSEA)')
plt.ylabel('')
plt.legend(loc=0, title='Type')
plt.xticks(np.arange(-1, 1.5, .5))
plt.gcf().set_size_inches(1, 3)
plt.title('GO terms')
plt.savefig('./reports/protein_enrichment_goterms.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Enrichment pathways
# Import
e_p_path = read_csv('./files/gsea_proteomics_metabolic_pathways.csv', index_col=2)
e_p_path = e_p_path[e_p_path['length'] >= 5]

e_f_path = read_csv('./files/gsea_fluxomics_metabolic_pathways.csv', index_col=2)
e_f_path = e_f_path[e_f_path['length'] >= 5]

# Overlap
ov_path = set(e_p_path.index).intersection(e_f_path.index)

# FDR
e_p_path = e_p_path.ix[ov_path]
e_p_path['fdr'] = multipletests(e_p_path['pvalue'], method='fdr_bh')[1]
e_p_path['abs_score'] = e_p_path['escore'].abs()
print e_p_path.sort('escore')

e_f_path = e_f_path.ix[ov_path]
e_f_path['fdr'] = multipletests(e_f_path['pvalue'], method='fdr_bh')[1]
e_f_path['abs_score'] = e_f_path['escore'].abs()
print e_f_path.sort('escore')

# Plot
plot_df = DataFrame([
    {
        'flux': e_f_path.ix[p, 'escore'],
        'protein': e_p_path.ix[p, 'escore'],
        'pathway': p
    } for p in ov_path
]).dropna().set_index('pathway')
plot_df = plot_df[plot_df['flux'].abs() > 1e-5]
print plot_df.sort('flux')

cor, pval = pearsonr(plot_df['flux'], plot_df['protein'])
print cor, pval

p_highlight = [
    'Citric acid cycle',
    'Glutamate metabolism',
    'Pyruvate metabolism',
    'Glycine, serine, alanine and threonine metabolism',
    'Pentose phosphate pathway',
    'beta-Alanine metabolism',
    'Glutathione metabolism'
]
pal = dict(zip(*(p_highlight + ['Others'], sns.color_palette('Set2', n_colors=8).as_hex() + ['#dfdfdf'])))

# Plot
sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
# sns.regplot('flux', 'protein', plot_df, fit_reg=False, color='#34495e', line_kws={'lw': .3})
for p in plot_df.index:
    if p in p_highlight:
        sns.regplot('flux', 'protein', plot_df.ix[[p]], fit_reg=False, label=p, color=pal[p], scatter_kws={'s': 30, 'lw': 1.})
    else:
        sns.regplot('flux', 'protein', plot_df.ix[[p]], fit_reg=False, color=pal['Others'], scatter_kws={'s': 10, 'lw': .0, 'edgecolors': 'none'})

plt.axhline(0, ls='-', lw=.3, c='gray')
plt.axvline(0, ls='-', lw=.3, c='gray')
plt.xlabel('Flux enrichment score')
plt.ylabel('Protein enrichment score')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Pearson r: %.2f, p-value: %.2e' % (cor, pval))
plt.gcf().set_size_inches(2, 2)
plt.savefig('./reports/protein_flux_scatter.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
