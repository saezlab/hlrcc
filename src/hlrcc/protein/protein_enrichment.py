#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import matplotlib.pyplot as plt
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv
from hlrcc.utils import read_gmt, get_complexes_dict, get_complexes_name


# Import proteomics
proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')

# Import GO terms
signatures = {
    'BP': read_gmt('./files/c5.bp.v5.1.symbols.gmt'),
    'CC': read_gmt('./files/c5.cc.v5.1.symbols.gmt'),
    'MF': read_gmt('./files/c5.mf.v5.1.symbols.gmt')
}

# -- Enrichment analysis
dataset = proteomics.to_dict()

# GO terms enrichment
df_enrichment = [(t, sig, len(db[sig].intersection(dataset)), gsea(dataset, db[sig], 1000)) for t, db in signatures.items() for sig in db]
df_enrichment = DataFrame([{'type': t, 'signature': s, 'length': l, 'escore': es, 'pvalue': pval} for t, s, l, (es, pval) in df_enrichment]).dropna()
df_enrichment.sort(['escore']).to_csv('./files/proteomics_tmt_go_term.csv', index=False)
# df_enrichment = read_csv('./files/proteomics_tmt_go_term.csv')
print df_enrichment[df_enrichment['length'] > 5].sort(['escore'])

# Protein complexes enrichment
corum, corum_n = get_complexes_dict(), get_complexes_name()

df_enrichment_complexes = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, 1000)) for c, sig in corum.items()]
df_enrichment_complexes = DataFrame([{'complex': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment_complexes]).dropna()
df_enrichment_complexes['name'] = [corum_n[i] for i in df_enrichment_complexes['complex']]
df_enrichment_complexes.sort(['escore']).to_csv('./files/proteomics_tmt_go_term_corum.csv', index=False)
# df_enrichment_complexes = read_csv('./files/proteomics_tmt_go_term_corum.csv')
print df_enrichment_complexes[df_enrichment_complexes['length'] > 4].sort(['escore'])
