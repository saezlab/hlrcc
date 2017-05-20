#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import re
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.stats.misc import zscore
from scipy.stats.stats import ttest_ind, ttest_1samp
from pandas import DataFrame, read_csv, concat
from statsmodels.stats.multitest import multipletests

# -- ID maps
umap = pickle.load(open('./files/uniprot_to_genename.pickle', 'rb'))


# -- Proteomics
# - Import samplesheet
samplesheet = read_csv('./data/proteomics_samplesheet.txt', sep='\t', index_col=0)
samplesheet = samplesheet.loc[np.bitwise_and(samplesheet['organism'] == 'human', samplesheet['type'] == 'tp')]

ko = samplesheet[samplesheet['condition'] == 'fh_ko'].index
wt = samplesheet[samplesheet['condition'] == 'fh_wt'].index

# - Import and process phospho
info_columns = ['peptide', 'uniprot']

proteomics = read_csv('./data/uok262_proteomics.txt', sep='\t').dropna(subset=info_columns)
proteomics = proteomics[np.concatenate((info_columns, samplesheet.index))].replace(0.0, np.NaN)

# - Remove peptides with 1 or less measurements per condition
proteomics = proteomics[np.bitwise_and(proteomics[ko].count(1) > 1, proteomics[wt].count(1) > 1)]

# - Considering proteotypic peptides
proteomics = proteomics[[(len(i.split('; ')) == 2) and (i.split('; ')[0] != '') for i in proteomics['uniprot']]]

# - Log 2 transform
proteomics[samplesheet.index] = np.log2(proteomics[samplesheet.index])

# - Scale samples
proteomics[samplesheet.index] = zscore(proteomics[samplesheet.index])


# Clustermap
cod_cmap = {'UOK262': u'#34495e', 'UOK262pFH': u'#919daa'}
cod_map = {'fh_ko': 'UOK262', 'fh_wt': 'UOK262pFH'}

plot_df = proteomics[samplesheet.index].corr(method='spearman').astype(float)
plot_df.index = [cod_map[samplesheet.loc[i, 'condition']] for i in plot_df.index]
plot_df.columns = [cod_map[samplesheet.loc[i, 'condition']] for i in plot_df]

row_color = DataFrame({'Condition': {i: cod_cmap[i] for i in plot_df.index}})
col_color = DataFrame({'Condition': {i: cod_cmap[i] for i in plot_df}})

cmap = sns.light_palette('#e65245', as_cmap=True)
sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, cmap=cmap, lw=.3, annot=True, fmt='.2f', row_colors=row_color, col_colors=col_color, figsize=(3, 3))
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.savefig('./reports/clustermap_proteomics.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'


# - Gene symbol map
proteomics = proteomics[[i.split(';')[0] in umap for i in proteomics['uniprot']]]
proteomics['genesymbol'] = [umap[i.split(';')[0]] for i in proteomics['uniprot']]

# - Log fold-change
proteomics = proteomics.groupby('genesymbol').mean()

# - Differential protein abundance
de_proteomics = {}
for i in proteomics.index:
    t, p = ttest_ind(proteomics.ix[i, ko], proteomics.ix[i, wt])
    de_proteomics[i] = {
        'fc': proteomics.ix[i, ko].mean() - proteomics.ix[i, wt].mean(),
        't': t,
        'pval': p
    }
de_proteomics = DataFrame(de_proteomics).T.dropna()

# - FDR correction
de_proteomics['fdr'] = multipletests(de_proteomics['pval'], method='fdr_bh')[1]

# - Export protein level proteomics
de_proteomics.to_csv('./data/uok262_proteomics_labelfree_processed_fc.csv')
print de_proteomics.sort_values('fdr')


# Volcano
plot_df = de_proteomics.copy()
plot_df['signif'] = ['*' if (i < 0.05) and (abs(f) > .5) else '-' for i, f in plot_df[['fdr', 'fc']].values]
plot_df['log_pval'] = -np.log10(plot_df['pval'])

sns.set(
    style='ticks', context='paper', font_scale=0.75,
    rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'}
)
pal = dict(zip(*(['*', '-'], sns.light_palette('#34495e', 3, reverse=True).as_hex()[:-1])))

g = sns.lmplot(
    x='fc', y='log_pval', data=plot_df, hue='signif', fit_reg=False, palette=pal, legend=False,
    scatter_kws={'edgecolor': 'w', 'lw': .3, 's': 20}
)
g.axes[0, 0].set_ylim(0,)
g.despine(right=False, top=False)

# Add FDR threshold lines
plt.text(plt.xlim()[0]*.98, -np.log10(plot_df.loc[plot_df['fdr'] < 0.05, 'pval'].max()) * 1.01, 'FDR 5%', ha='left', color=pal['*'], alpha=0.65, fontsize=5)
plt.axhline(-np.log10(plot_df.loc[plot_df['fdr'] < 0.05, 'pval'].max()), c=pal['*'], ls='--', lw=.3, alpha=.7)

plt.axvline(-.5, c=pal['*'], ls='--', lw=.3, alpha=.7)
plt.axvline(.5, c=pal['*'], ls='--', lw=.3, alpha=.7)

# Add axis lines
plt.axvline(0, c='#95a5a6', lw=.3, ls='-', alpha=.3)

# Add axis labels and title
plt.title('Proteomics fold-change\n(UOK262 - UOK262pFH)', fontsize=8, fontname='sans-serif')
plt.xlabel('Fold-change (log2)', fontsize=8, fontname='sans-serif')
plt.ylabel('p-value (-log10)', fontsize=8, fontname='sans-serif')

# Add text to highlighted genes
for g in ['FH', 'VIM', 'GAPDH']:
    plt.scatter([plot_df.ix[g, 'fc']], [plot_df.ix[g, 'log_pval']], color='#e74c3c', alpha=.8, s=15)
    plt.text((plot_df.ix[g, 'fc'] * 1.01), (plot_df.ix[g, 'log_pval'] * 1.01), g, ha='left', alpha=0.75, fontsize=8, fontname='sans-serif')

# Save plot
plt.gcf().set_size_inches(3., 5., forward=True)
plt.savefig('./reports/volcano_proteomics.pdf', bbox_inches='tight')
plt.close('all')


# -- Phosphoproteomics
samplesheet = read_csv('./data/proteomics_samplesheet.txt', sep='\t', index_col=0)
samplesheet = samplesheet.loc[np.bitwise_and(samplesheet['organism'] == 'human', samplesheet['type'] == 'pp')]

ko = samplesheet[samplesheet['condition'] == 'fh_ko'].index
wt = samplesheet[samplesheet['condition'] == 'fh_wt'].index

# - Import and process phospho
info_columns = ['peptide', 'site', 'uniprot']

phosphoproteomics = read_csv('./data/uok262_phosphoproteomics.txt', sep='\t').dropna(subset=info_columns)
phosphoproteomics = phosphoproteomics[info_columns + list(samplesheet.index)].replace(0.0, np.NaN)

# - Remove peptides with 1 or less measurements per condition
phosphoproteomics = phosphoproteomics[(phosphoproteomics[ko].count(1) > 1) & (phosphoproteomics[wt].count(1) > 1)]

# - Log 2 transform
phosphoproteomics[samplesheet.index] = np.log2(phosphoproteomics[samplesheet.index])

# - Scale samples
phosphoproteomics[samplesheet.index] = zscore(phosphoproteomics[samplesheet.index])


# Clustermap
cod_cmap = {'UOK262': u'#34495e', 'UOK262pFH': u'#919daa'}
cod_map = {'fh_ko': 'UOK262', 'fh_wt': 'UOK262pFH'}

plot_df = phosphoproteomics[samplesheet.index].corr(method='spearman').astype(float)
plot_df.index = [cod_map[samplesheet.loc[i, 'condition']] for i in plot_df.index]
plot_df.columns = [cod_map[samplesheet.loc[i, 'condition']] for i in plot_df]

row_color = DataFrame({'Condition': {i: cod_cmap[i] for i in plot_df.index}})
col_color = DataFrame({'Condition': {i: cod_cmap[i] for i in plot_df}})

cmap = sns.light_palette('#e65245', as_cmap=True)
sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, cmap=cmap, lw=.3, annot=True, fmt='.2f', row_colors=row_color, col_colors=col_color, figsize=(4, 4))
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.savefig('./reports/clustermap_phosphoproteomics.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'


# - Considering proteotypic peptides
phosphoproteomics = phosphoproteomics[[len(i.split('; ')) == 2 for i in phosphoproteomics['uniprot']]]

# - Consider single phosphorylated peptides
phosphoproteomics = phosphoproteomics[[len(i.split('+')) == 1 for i in phosphoproteomics['site']]]

# - Map uniprot to gene symbol
phosphoproteomics = phosphoproteomics[[i.split(';')[0] in umap for i in phosphoproteomics['uniprot']]]
phosphoproteomics['genesymbol'] = [umap[i.split(';')[0]] for i in phosphoproteomics['uniprot']]

# - Create p-site IDs
phosphoproteomics['psite'] = ['%s_%s' % (g, re.findall('\(([A-Z][0-9]*)\)', p)[0]) for g, p in phosphoproteomics[['genesymbol', 'site']].values]

# - Log fold-change
phosphoproteomics = phosphoproteomics.groupby('psite').mean()

# - Differential protein abundance
de_phosphoproteomics = {}
for i in phosphoproteomics.index:
    t, p = ttest_ind(phosphoproteomics.ix[i, ko], phosphoproteomics.ix[i, wt])
    de_phosphoproteomics[i] = {
        'fc': phosphoproteomics.ix[i, ko].mean() - phosphoproteomics.ix[i, wt].mean(),
        't': t,
        'pval': p
    }
de_phosphoproteomics = DataFrame(de_phosphoproteomics).T.dropna()

# - FDR correction
de_phosphoproteomics['fdr'] = multipletests(de_phosphoproteomics['pval'], method='fdr_bh')[1]

# - Export p-site level phosphoproteomics
de_phosphoproteomics.to_csv('./data/uok262_phosphoproteomics_labelfree_processed_fc.csv')
print de_phosphoproteomics.sort_values('fdr')


# - Volcano
plot_df = de_phosphoproteomics.copy()
plot_df['signif'] = ['*' if (i < 0.05) and (abs(f) > .5) else '-' for i, f in plot_df[['fdr', 'fc']].values]
plot_df['log_pval'] = -np.log10(plot_df['pval'])

sns.set(
    style='ticks', context='paper', font_scale=0.75,
    rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'}
)
pal = dict(zip(*(['*', '-'], sns.light_palette('#34495e', 3, reverse=True).as_hex()[:-1])))

g = sns.lmplot(
    x='fc', y='log_pval', data=plot_df, hue='signif', fit_reg=False, palette=pal, legend=False,
    scatter_kws={'edgecolor': 'w', 'lw': .3, 's': 20}
)
g.axes[0, 0].set_ylim(0,)
g.despine(right=False, top=False)

# Add FDR threshold lines
plt.text(plt.xlim()[0]*.98, -np.log10(plot_df.loc[plot_df['fdr'] < 0.05, 'pval'].max()) * 1.01, 'FDR 5%', ha='left', color=pal['*'], alpha=0.65, fontsize=5)
plt.axhline(-np.log10(plot_df.loc[plot_df['fdr'] < 0.05, 'pval'].max()), c=pal['*'], ls='--', lw=.3, alpha=.7)

plt.axvline(-.5, c=pal['*'], ls='--', lw=.3, alpha=.7)
plt.axvline(.5, c=pal['*'], ls='--', lw=.3, alpha=.7)

# Add axis lines
plt.axvline(0, c='#95a5a6', lw=.3, ls='-', alpha=.3)

# Add axis labels and title
plt.title('Proteomics fold-change\n(UOK262 - UOK262pFH)', fontsize=8, fontname='sans-serif')
plt.xlabel('Fold-change (log2)', fontsize=8, fontname='sans-serif')
plt.ylabel('p-value (-log10)', fontsize=8, fontname='sans-serif')

# Add text to highlighted genes
for p in ['VIM', 'GAPDH', 'PDHA1']:
    for g in [i for i in plot_df.index if p in i]:
        plt.scatter([plot_df.ix[g, 'fc']], [plot_df.ix[g, 'log_pval']], color='#e74c3c', alpha=.8, s=15)
        plt.text((plot_df.ix[g, 'fc'] * 1.01), (plot_df.ix[g, 'log_pval'] * 1.01), g, ha='left', alpha=0.75, fontsize=8, fontname='sans-serif')

# Save plot
plt.gcf().set_size_inches(3., 5., forward=True)
plt.savefig('./reports/volcano_phosphoproteomics.pdf', bbox_inches='tight')
plt.close('all')
