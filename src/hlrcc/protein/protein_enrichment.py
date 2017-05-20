#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import matplotlib.pyplot as plt
from framed import load_cbmodel
from hlrcc.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv, concat
from hlrcc.utils import read_gmt, get_complexes_dict, get_complexes_name


# -- Import data-sets
proteomics = read_csv('./data/uok262_proteomics_labelfree_processed_fc.csv', index_col=0)['fc']
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')


# -- Import gene signatures
# GO terms
signatures = {
    'BP': read_gmt('./files/c5.bp.v5.1.symbols.gmt'),
    'CC': read_gmt('./files/c5.cc.v5.1.symbols.gmt'),
    'MF': read_gmt('./files/c5.mf.v5.1.symbols.gmt')
}

# sig = 'SERINE_TYPE_ENDOPEPTIDASE_INHIBITOR_ACTIVITY'
# gsea(dataset, signatures['MF'][sig], 1000, './reports/gsea_MF_plot.pdf', sig.lower().replace('_', ' ').capitalize())
#
# sig = 'MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_I'
# gsea(dataset, signatures['CC'][sig], 1000, './reports/gsea_CC_plot.pdf', sig.lower().replace('_', ' ').capitalize())

# Corum
corum, corum_n = get_complexes_dict(), get_complexes_name()

# Metabolic pathways
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])

model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')

r_genes = {r: {gmap.ix[g] for g in model.reactions[r].gpr.get_genes() if g in gmap.index} for r in model.reactions if model.reactions[r].gpr}

p_reactions = DataFrame([{'r': r, 'p': model.reactions[r].metadata['SUBSYSTEM']} for r in model.reactions if 'SUBSYSTEM' in model.reactions[r].metadata and model.reactions[r].metadata['SUBSYSTEM'] != ''])
p_reactions = p_reactions.groupby('p')['r'].agg(lambda x: set(x)).to_dict()
p_reactions_genes = {p: {g for r in p_reactions[p] if r in r_genes for g in r_genes[r]} for p in p_reactions}


# -- Enrichment analysis
nrand, dataset = 1000, proteomics.to_dict()

# GO terms enrichment
df_enrichment = [(t, sig, len(db[sig].intersection(dataset)), gsea(dataset, db[sig], nrand)) for t, db in signatures.items() for sig in db]
df_enrichment = DataFrame([{'type': t, 'signature': s, 'length': l, 'escore': es, 'pvalue': pval} for t, s, l, (es, pval) in df_enrichment]).dropna()
df_enrichment.sort(['escore']).to_csv('./files/gsea_proteomics_goterms.csv', index=False)
print df_enrichment[df_enrichment['length'] >= 3].sort(['escore'])

# Protein complexes enrichment
df_enrichment_complexes = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, nrand)) for c, sig in corum.items()]
df_enrichment_complexes = DataFrame([{'complex': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment_complexes]).dropna()
df_enrichment_complexes['name'] = [corum_n[i] for i in df_enrichment_complexes['complex']]
df_enrichment_complexes.sort(['escore']).to_csv('./files/gsea_proteomics_corum.csv', index=False)
print df_enrichment_complexes[df_enrichment_complexes['length'] >= 5].sort(['escore'])

# Metabolic pathways enrichment
df_enrichment_pathways = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, nrand)) for c, sig in p_reactions_genes.items()]
df_enrichment_pathways = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment_pathways]).dropna()
df_enrichment_pathways.sort(['escore']).to_csv('./files/gsea_proteomics_metabolic_pathways.csv', index=False)
print df_enrichment_pathways[df_enrichment_pathways['length'] >= 3].sort(['escore'])


# -- Import fluxomics
conditions = ['KO', 'WT']
conditions_map = {'UOK262': 'KO', 'UOK262pFH': 'WT'}

fluxes = read_csv('./data/pfba_atp.csv', index_col=0).replace(np.nan, 0)
fluxes['delta'] = fluxes['UOK262'] - fluxes['UOK262pFH']
fluxes.columns = [conditions_map[c] if c in conditions_map else c for c in fluxes]

# Enrichment
nrand, dataset = 1000, fluxes['delta'].to_dict()

# Metabolic pathways enrichment
df_flux_enrichment_pathways = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, nrand)) for c, sig in p_reactions.items()]
df_flux_enrichment_pathways = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_flux_enrichment_pathways]).dropna()
df_flux_enrichment_pathways.sort(['escore']).to_csv('./files/gsea_fluxomics_metabolic_pathways.csv', index=False)
print df_flux_enrichment_pathways[df_flux_enrichment_pathways['length'] >= 3].sort(['escore'])
