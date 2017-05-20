#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import pickle
from pandas import read_csv

# - Uniprot ID mapping
umap = read_csv('./data/HUMAN_9606_idmapping.txt', sep='\t')
umap = umap[umap['database'] == 'Gene_Name']
umap = umap.groupby('uniprot')['id'].agg(lambda x: set(x)).to_dict()
umap = {k: list(umap[k])[0] for k in umap if len(umap[k]) == 1}

# - Store as a pickle object
with open('./files/uniprot_to_genename.pickle', 'wb') as handle:
    pickle.dump(umap, handle, protocol=pickle.HIGHEST_PROTOCOL)
