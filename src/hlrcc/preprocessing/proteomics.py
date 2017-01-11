#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv
from pandas.stats.misc import zscore, DataFrame
from hlrcc.utils.volcano import volcano
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.reader.sbml_reader import read_sbml_model


# -- Uniprot names
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)


# -- Import metabolic model
m_genes = read_sbml_model('./files/recon1.xml').get_genes()


# -- Import samplesheet
ss = read_csv('./data/proteomics_samplesheet.txt', sep='\t', index_col=0)
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'tp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index


# -- Process
# Import and process phospho
info_columns = ['peptide', 'uniprot']

tp_all = read_csv('./data/uok262_proteomics.txt', sep='\t').dropna(subset=info_columns)
tp = tp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', tp.shape

# Remove peptides with 1 or less measurements per condition
tp = tp[np.bitwise_and(tp[ss_ko].count(1) > 1, tp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', tp.shape

# Considering proteotypic peptides
tp = tp[[(len(i.split('; ')) == 2) and i.split('; ')[0] != '' for i in tp['uniprot']]]
print '[INFO] Considering proteotypic peptides: ', tp.shape

# Create protein IDs
tp['uniprot'] = [i.split(';')[0] for i in tp['uniprot']]
tp = tp.groupby('uniprot', sort=False).median()

# Log 2 transform
tp[ss.index] = np.log2(tp[ss.index])

# Scale samples
tp[ss.index] = zscore(tp[ss.index])

# Export protein level proteomics
tp.columns = [ss.ix[i, 'condition'] for i in tp.columns]
tp.to_csv('./data/uok262_proteomics_processed.txt', sep='\t')
print '[INFO] Export protein level proteomics'


# -- Plot
# Heatmap
cmap = sns.light_palette((210, 90, 60), input='husl', as_cmap=True)

# df = DataFrame(np.triu(tp.corr(method='pearson'), 1), columns=tp.columns, index=tp.columns).replace(0, np.nan).unstack().reset_index().dropna()
# df[df['level_0'] == df['level_1']][0].mean()
# df[df['level_0'] != df['level_1']][0].mean()

sns.clustermap(tp.corr(method='pearson'), annot=True, cmap=cmap)
plt.savefig('./reports/proteomics_replicates_clustermap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'

# Volcano
tp_fc = read_csv('./data/uok262_proteomics_logfc.txt', sep='\t')
tp_fc['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else '' for i in tp_fc.index]

genes_highlight = ['VIM', 'PDHA1', 'GAPDH', 'FH']

volcano(
    './reports/proteomics_logfc_volcano.pdf',
    tp_fc,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Proteomics (KO vs WT)',
    genes_highlight,
    'name'
)
plt.close('all')
print '[INFO] Plot done!'
