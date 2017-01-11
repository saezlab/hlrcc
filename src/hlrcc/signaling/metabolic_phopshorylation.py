#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from pymist.reader.sbml_reader import read_sbml_model
from pandas import DataFrame, read_csv
from statsmodels.stats.multitest import multipletests
from pymist.utils.map_peptide_sequence import read_uniprot_genename


def cohensd(c0, c1):
    # return (np.mean(c0) - np.mean(c1)) / (np.sqrt((np.std(c0) ** 2 + np.std(c1) ** 2) / 2))
    return np.mean(c0) - np.mean(c1)


# -- Imports
# Uniprot id mapping
human_uniprot = read_uniprot_genename()

# Metabolic model
m_model = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml')
m_genes = m_model.get_genes()

# Phosphoproteomics fold-change
phospho_fc = read_csv('./data/uok262_phosphoproteomics_logfc.txt', sep='\t')
phospho_fc['psite'] = ['_'.join(i.split('_')[:2]) for i in phospho_fc.index]
phospho_fc = phospho_fc[[i.split('_')[0] in human_uniprot for i in phospho_fc['psite']]]
phospho_fc['psite'] = ['%s_%s' % (human_uniprot[i.split('_')[0]][0], i.split('_')[1]) for i in phospho_fc['psite']]

# Metabolism sampling
ko_sampling, wt_sampling = [read_csv('./data/%s_sampling.txt' % c, sep='\t', index_col=0) for c in ['UOK262', 'UOK262pFH']]

metab_fc = [(r, cohensd(ko_sampling[r], wt_sampling[r]), wilcoxon(ko_sampling[r], wt_sampling[r])[1]) for r in ko_sampling]
metab_fc = DataFrame(metab_fc, columns=['reaction', 'cohensd', 'wilcoxon'])
metab_fc = metab_fc[metab_fc['wilcoxon'] != 0]
metab_fc['adj.P.Val'] = multipletests(metab_fc['wilcoxon'], method='fdr_bh')[1]
metab_fc = metab_fc.set_index('reaction')


# -- Overlap phosphorylation changes with metabolic changes
enzymes_phospho = phospho_fc[[i.split('_')[0] in m_genes for i in phospho_fc['psite']]]

enzymes_metabol = enzymes_phospho[['psite', 'logFC', 'adj.P.Val']]
enzymes_metabol['enzyme'] = [i.split('_')[0] for i in enzymes_metabol['psite']]
enzymes_metabol = [(enz, psite, logfc, pval, r, metab_fc.ix[r, 'cohensd'], metab_fc.ix[r, 'adj.P.Val']) for psite, logfc, pval, enz in enzymes_metabol.values for r in m_model.get_reactions_by_genes([enz])[enz] if r in metab_fc.index]
enzymes_metabol = DataFrame(enzymes_metabol, columns=['enzyme', 'psite', 'psite_logfc', 'psite_fdr', 'reaction', 'reaction_cohensd', 'reaction_fdr'])
enzymes_metabol = enzymes_metabol[(enzymes_metabol['psite_fdr'] < 0.05) & (enzymes_metabol['reaction_fdr'] < 0.05)]
enzymes_metabol['reaction_cohensd_abs'] = enzymes_metabol['reaction_cohensd'].abs()
enzymes_metabol = enzymes_metabol.sort(['reaction_cohensd_abs'], ascending=False)
print enzymes_metabol


# -- Plot
# Metabolic enzymes phosphorylation
pal = sns.light_palette('#34495e', len(enzymes_metabol) + 1, reverse=True)[:-1]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.lmplot('psite_logfc', 'reaction_cohensd', enzymes_metabol, 'psite', fit_reg=False, palette=pal, scatter_kws={'s': 50, 'edgecolor': 'w', 'linewidth': .5})
plt.axhline(0, ls='-', lw=.3, c='gray')
plt.axvline(0, ls='-', lw=.3, c='gray')
plt.xlabel('Phosphorylation-site (log2 FC)')
plt.ylabel('Flux rate (mean difference)')
plt.title('UOK262 (KO vs WT)')
plt.gcf().set_size_inches(5, 4)
plt.savefig('./reports/psites_reactions_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'

# Metabolic sampling
reactions = list(enzymes_metabol['reaction'])

r_enz = enzymes_metabol[[i in reactions for i in enzymes_metabol['reaction']]][['reaction', 'enzyme']].set_index('reaction').to_dict()['enzyme']

ko_samples = ko_sampling[reactions].unstack().reset_index()
ko_samples['condition'] = 'UOK262'

wt_samples = wt_sampling[reactions].unstack().reset_index()
wt_samples['condition'] = 'UOK262pFH'

plot_df = ko_samples.append(wt_samples)
plot_df.columns = ['reaction', 'sample', 'value', 'condition']
plot_df['enzyme'] = [r_enz[i] for i in plot_df['reaction']]
plot_df['reaction'] = [i[2:] for i in plot_df['reaction']]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, col='enzyme', col_wrap=4, sharey=False, sharex=False, aspect=1.7, size=1.3, legend_out=True)
g.map(sns.violinplot, 'value', 'reaction', 'condition', orient='h', palette=sns.light_palette('#34495e', 3)[1:], split=True, inner='quart', cut=0)
sns.despine(trim=True)
g.set_titles('{col_name}')
g.set_ylabels('')
g.set_xlabels('Flux rate (umol/ugDW/h)')
g.add_legend()
plt.savefig('./reports/sampling_boxplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
