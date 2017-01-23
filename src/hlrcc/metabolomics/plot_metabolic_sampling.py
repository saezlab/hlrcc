#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from framed import load_cbmodel
from pandas import DataFrame, Series, read_csv
from statsmodels.stats.multitest import multipletests
from more_itertools import unique_everseen

# -- Imports
# Samples
ko_sampling, wt_sampling = [read_csv('./data/%s_sampling.txt' % c, sep='\t') for c in ['UOK262', 'UOK262pFH']]

# Model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Phosphoproteomics
umap = read_csv('./files/protein-coding_gene.txt', sep='\t').dropna(subset=['uniprot_ids'])
umap = umap.groupby('uniprot_ids')['symbol'].agg(lambda x: ';'.join([g for i in x for g in i.split('|')]))

phospho_fc = read_csv('./data/uok262_phosphoproteomics_logfc.txt', sep='\t')
phospho_fc['psite'] = ['_'.join(i.split('_')[:2]) for i in phospho_fc.index]
phospho_fc = phospho_fc[[i.split('_')[0] in umap.index for i in phospho_fc['psite'] if i.split('_')[0]]]
phospho_fc['gene'] = [umap.ix[i.split('_')[0]] for i in phospho_fc['psite']]
print phospho_fc[phospho_fc['adj.P.Val'] < .05].sort('logFC')


# -- Reactions map
gmap = read_csv('./files/protein-coding_gene.txt', sep='\t')
gmap['hgnc_id'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.set_index('hgnc_id')['symbol']

r_genes = {r: {gmap.ix[g] for g in model.reactions[r].gpr.get_genes() if g in gmap.index} for r in model.reactions if model.reactions[r].gpr}


# --
metab_fc = DataFrame(
    {r: {'UOK262': ko_sampling[r].mean(), 'UOK262pFH': wt_sampling[r].mean(), 'wilcoxon': wilcoxon(ko_sampling[r], wt_sampling[r])[1]} for r in ko_sampling}
).T
metab_fc['diff'] = metab_fc['UOK262'] - metab_fc['UOK262pFH']
metab_fc['fdr'] = multipletests(metab_fc['wilcoxon'], method='bonferroni')[1]
metab_fc['genes'] = [';'.join(r_genes[r]) if r in r_genes else np.nan for r in metab_fc.index]
metab_fc['pathway'] = [model.reactions[r].metadata['SUBSYSTEM'] if 'SUBSYSTEM' in model.reactions[r].metadata else np.nan for r in metab_fc.index]
print metab_fc[(metab_fc['fdr'] < .05)].dropna().sort('diff')


# --
pgenes = set(phospho_fc.ix[(phospho_fc['adj.P.Val'] < .05) & (phospho_fc['logFC'].abs() > .5), 'gene'])

link = metab_fc[(metab_fc['fdr'] < .05) & (metab_fc['diff'].abs() > .5)].copy().dropna(subset=['genes']).reset_index()
link = link[[len(set(i.split(';')).intersection(pgenes)) > 0 for i in link['genes']]]
link = link[link['genes'] != 'CMPK1']

link = DataFrame([
    {'flux': d, 'reaction': r, 'complex': gs, 'gene': g, 'psite': '_'.join([g, psite.split('_')[1]]), 'phospho': logfc}
    for d, gs, r in link[['diff', 'genes', 'index']].values for g in gs.split(';') for psite, logfc in phospho_fc.ix[phospho_fc['gene'] == g, ['psite', 'logFC']].values
]).sort('phospho')
print link

# Plot
order = list(unique_everseen(link['psite']))
pal = dict(zip(*(order, sns.light_palette('#34495e', len(order)+1, reverse=True).as_hex()[:-1])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.lmplot('phospho', 'flux', link, 'psite', fit_reg=False, palette=pal, scatter_kws={'s': 50, 'edgecolor': 'w', 'linewidth': .5}, size=3)
plt.axhline(0, ls='-', lw=.3, c='gray')
plt.axvline(0, ls='-', lw=.3, c='gray')
plt.xlabel('Phosphorylation-site (log2 FC)')
plt.ylabel('Flux (mean difference)')
plt.title('UOK262 (KO vs WT)')
plt.savefig('./reports/psites_reactions_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'

# -- Plot
reactions = ['R_OIVD2m', 'R_OIVD3m', 'R_GAPD']

ko_samples = ko_sampling[reactions].unstack().reset_index()
ko_samples.columns = ['reaction', 'index', 'flux']
ko_samples['condition'] = 'UOK262'

wt_samples = wt_sampling[reactions].unstack().reset_index()
wt_samples.columns = ['reaction', 'index', 'flux']
wt_samples['condition'] = 'UOK262pFH'

plot_df = ko_samples.append(wt_samples)

pal = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#34495e', 3)[1:])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, row='reaction', sharey=False, sharex=False, aspect=.75, size=1.5, legend_out=True)
g.map(sns.violinplot, 'flux', 'condition', orient='h', palette=pal, split=True, inner='quart', cut=0, lw=.3)
sns.despine(trim=True)
# g.set_titles('{col_name}')
g.set_ylabels('')
g.set_xlabels('Flux rate (mmol/gDW/h)')
g.add_legend()
plt.savefig('./reports/sampling_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
