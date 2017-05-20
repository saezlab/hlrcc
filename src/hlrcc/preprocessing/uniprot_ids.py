#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import pickle
from pandas import read_csv

# Convert UniProt ids to Gene Symbols
umap = read_csv('./data/HUMAN_9606_idmapping.txt', sep='\t')

# Gene symbol only
umap = umap[umap['database'] == 'Gene_Name']

# Build dict
umap = umap.groupby('uniprot')['id'].agg(lambda x: list(x)).to_dict()

# Exclude ambigous mapping
umap = {k: umap[k][0] for k in umap if len(umap[k]) == 1}

# Store as a pickle object
with open('./data/uniprot_to_genename.pickle', 'wb') as handle:
    pickle.dump(umap, handle, protocol=pickle.HIGHEST_PROTOCOL)
