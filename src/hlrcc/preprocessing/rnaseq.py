#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import matplotlib.pyplot as plt
from pandas import Series, read_csv

# -- Gene-name map
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap = gmap.groupby('ensembl_gene_id')['symbol'].agg(lambda x: list(x))


# -- Import RNA-seq
rnaseq = read_csv('./data/UOK262_rnaseq.csv', index_col=0)
rnaseq = rnaseq[[g in gmap and len(gmap[g]) == 1 for g in rnaseq.index]]
rnaseq['symbol'] = [gmap[g][0] for g in rnaseq.index]

rnaseq_fc = rnaseq.set_index('symbol')['logFC'].copy()
rnaseq_fc.to_csv('./data/UOK262_rnaseq_preprocessed.csv')
print rnaseq_fc.sort_values()
