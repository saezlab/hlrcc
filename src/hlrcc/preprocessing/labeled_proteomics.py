#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import re
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc.utils import read_fasta
from scipy.stats.stats import ttest_1samp
from pandas import DataFrame, read_csv, concat
from statsmodels.stats.multitest import multipletests


def map_sites(acc, pep, ptms):
    if acc in human_uniprot and len(re.findall(pep, human_uniprot[acc])) == 1:
        pep_start = human_uniprot[acc].find(pep)
        return [ptm[0] + str(pep_start + int(re.findall('([0-9]+)\(', ptm)[0])) for ptm in ptms]

    else:
        return []


# -- Import uniprot sequences
human_uniprot = read_fasta()


# -- ID maps
umap = pickle.load(open('./files/uniprot_to_genename.pickle', 'rb'))


# -- Process TMT proteomics
proteomics = read_csv('./data/uok262_proteomics_tmt.csv', index_col=0)[['128_C/127_N', '128_N/126', '129_N/127_C']].dropna()

# Log2
proteomics = np.log2(proteomics)

# Map gene name
proteomics = proteomics[[i in umap for i in proteomics.index]]
proteomics['gene'] = [umap[i] for i in proteomics.index]

# Average peptides matching to same protein
proteomics = proteomics.groupby('gene').median()

# Differential protein abundance
de_proteomics = {}
for g in proteomics.index:
    t, p = ttest_1samp(proteomics.ix[g], 0)
    de_proteomics[g] = {'fc': proteomics.ix[g].mean(), 't': t, 'pval': p}
de_proteomics = concat([proteomics, DataFrame(de_proteomics).T], axis=1)

# FDR correction
de_proteomics['fdr'] = multipletests(de_proteomics['pval'], method='fdr_bh')[1]

# Export
proteomics.to_csv('./data/uok262_proteomics_tmt_preprocessed.csv')
print '[INFO] TMT proteomics exported'


# -- Process TMT phosphoproteomics
phosphoproteomics = read_csv('./data/uok262_phosphoproteomics_tmt.csv').dropna(subset=['Sequence', 'Protein Group Accessions', 'Modifications'])

# Consider only phospho PTMs
phosphoproteomics['Modifications'] = ['; '.join([ptm for ptm in m.split('; ') if 'Phospho' in ptm]) for m in phosphoproteomics['Modifications']]

# Upper case sequences
phosphoproteomics['Sequence'] = [s.upper() for s in phosphoproteomics['Sequence']]

# Map protein p-site
phosphoproteomics['sites'] = [map_sites(acc, pep, ptms.split('; ')) for pep, acc, ptms in phosphoproteomics[['Sequence', 'Protein Group Accessions', 'Modifications']].values]
phosphoproteomics = phosphoproteomics[[len(s) != 0 for s in phosphoproteomics['sites']]]

# Uniprot to gene symbol
phosphoproteomics = phosphoproteomics[[i in umap for i in phosphoproteomics['Protein Group Accessions']]]
phosphoproteomics['gene'] = [umap[i] for i in phosphoproteomics['Protein Group Accessions']]

# Create psite id
phosphoproteomics['sites'] = [acc + '_' + '_'.join(sites) for acc, sites in phosphoproteomics[['gene', 'sites']].values]

# Log2
phosphoproteomics.ix[:, ['128_C/127_N', '128_N/126', '129_N/127_C']] = np.log2(phosphoproteomics[['128_C/127_N', '128_N/126', '129_N/127_C']])

# Average p-sites within sample
phosphoproteomics = phosphoproteomics.groupby('sites')['128_C/127_N', '128_N/126', '129_N/127_C'].median().dropna()

# Differential phosphorylation
de_phosphoproteomics = {}
for g in phosphoproteomics.index:
    t, p = ttest_1samp(phosphoproteomics.ix[g], 0)
    de_phosphoproteomics[g] = {'fc': phosphoproteomics.ix[g].mean(), 't': t, 'pval': p}
de_phosphoproteomics = concat([phosphoproteomics, DataFrame(de_phosphoproteomics).T], axis=1)

# FDR correction
de_phosphoproteomics['fdr'] = multipletests(de_phosphoproteomics['pval'], method='fdr_bh')[1]

# Export
de_phosphoproteomics.to_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')
print '[INFO] TMT phosphoproteomics exported'
