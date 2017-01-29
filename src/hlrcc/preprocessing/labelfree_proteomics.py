#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import re
import pickle
import numpy as np
from pandas import read_csv
from pandas.stats.misc import zscore

# -- ID maps
with open('./files/uniprot_to_genename.pickle', 'rb') as handle:
    gmap = pickle.load(handle)

# -- Proteomics
# Import samplesheet
samplesheet = read_csv('./data/proteomics_samplesheet.txt', sep='\t', index_col=0)
samplesheet = samplesheet.loc[np.bitwise_and(samplesheet['organism'] == 'human', samplesheet['type'] == 'tp')]

ko = samplesheet[samplesheet['condition'] == 'fh_ko'].index
wt = samplesheet[samplesheet['condition'] == 'fh_wt'].index

# Import and process phospho
info_columns = ['peptide', 'uniprot']

proteomics = read_csv('./data/uok262_proteomics.txt', sep='\t').dropna(subset=info_columns)
proteomics = proteomics[np.concatenate((info_columns, samplesheet.index))].replace(0.0, np.NaN)

# Remove peptides with 1 or less measurements per condition
proteomics = proteomics[np.bitwise_and(proteomics[ko].count(1) > 1, proteomics[wt].count(1) > 1)]

# Considering proteotypic peptides
proteomics = proteomics[[(len(i.split('; ')) == 2) and (i.split('; ')[0] != '') for i in proteomics['uniprot']]]

# Log 2 transform
proteomics[samplesheet.index] = np.log2(proteomics[samplesheet.index])

# Scale samples
proteomics[samplesheet.index] = zscore(proteomics[samplesheet.index])

# Log fold-change
proteomics['logfc'] = proteomics[ko].mean(1) - proteomics[wt].mean(1)

# Gene symbol map
proteomics = proteomics[[i.split(';')[0] in gmap for i in proteomics['uniprot']]]
proteomics = proteomics[[len(gmap[i.split(';')[0]]) == 1 for i in proteomics['uniprot']]]
proteomics['genesymbol'] = [gmap[i.split(';')[0]][0] for i in proteomics['uniprot']]

# Fold-change
proteomics = proteomics.groupby('genesymbol')['logfc'].mean()

# Export protein level proteomics
proteomics.to_csv('./data/uok262_proteomics_labelfree_processed.csv')
print proteomics.sort_values()


# -- Phosphoproteomics
samplesheet = read_csv('./data/proteomics_samplesheet.txt', sep='\t', index_col=0)
samplesheet = samplesheet.loc[np.bitwise_and(samplesheet['organism'] == 'human', samplesheet['type'] == 'pp')]

ko = samplesheet[samplesheet['condition'] == 'fh_ko'].index
wt = samplesheet[samplesheet['condition'] == 'fh_wt'].index

# Import and process phospho
info_columns = ['peptide', 'site', 'uniprot']

phosphoproteomics = read_csv('./data/uok262_phosphoproteomics.txt', sep='\t').dropna(subset=info_columns)
phosphoproteomics = phosphoproteomics[info_columns + list(samplesheet.index)].replace(0.0, np.NaN)

# Remove peptides with 1 or less measurements per condition
phosphoproteomics = phosphoproteomics[(phosphoproteomics[ko].count(1) > 1) & (phosphoproteomics[wt].count(1) > 1)]

# Log 2 transform
phosphoproteomics[samplesheet.index] = np.log2(phosphoproteomics[samplesheet.index])

# Scale samples
phosphoproteomics[samplesheet.index] = zscore(phosphoproteomics[samplesheet.index])

# Considering proteotypic peptides
phosphoproteomics = phosphoproteomics[[len(i.split('; ')) == 2 for i in phosphoproteomics['uniprot']]]

# Consider single phosphorylated peptides
phosphoproteomics = phosphoproteomics[[len(i.split('+')) == 1 for i in phosphoproteomics['site']]]

# Map uniprot to gene symbol
phosphoproteomics = phosphoproteomics[[i.split(';')[0] in gmap for i in phosphoproteomics['uniprot']]]
phosphoproteomics = phosphoproteomics[[len(gmap[i.split(';')[0]]) == 1 for i in phosphoproteomics['uniprot']]]
phosphoproteomics['genesymbol'] = [gmap[i.split(';')[0]][0] for i in phosphoproteomics['uniprot']]

# Create p-site IDs
phosphoproteomics['psite'] = ['%s_%s' % (g, re.findall('\(([A-Z][0-9]*)\)', p)[0]) for g, p in phosphoproteomics[['genesymbol', 'site']].values]

# Log fold-change
phosphoproteomics['logfc'] = phosphoproteomics[ko].mean(1) - phosphoproteomics[wt].mean(1)

# Fold-change
phosphoproteomics = phosphoproteomics.groupby('psite')['logfc'].mean()

# Export p-site level phosphoproteomics
phosphoproteomics.to_csv('./data/uok262_phosphoproteomics_labelfree_processed.csv')
print phosphoproteomics.sort_values()
