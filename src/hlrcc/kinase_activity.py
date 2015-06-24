import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from hlrcc import wd
from pymist.enrichment.gsea import gsea
from scipy.stats.distributions import hypergeom
from bioservices import KEGG, KEGGParser
from pandas import DataFrame, Series, read_csv

# Import phospho data
pp = read_csv(wd + '/files/b1368p100_phospho_human_processed.tsv', sep='\t', index_col=0)

# Fold-change analysis KO - WT
ko, wt = [i for i in pp.columns if i.startswith('fh_ko')], [i for i in pp.columns if i.startswith('fh_wt')]
pp_fc = pp[ko].median(1) - pp[wt].median(1)

# ---- Calculate kinase activity scores
# Import network
ks_network = read_csv('/Users/emanuel/Projects/resources/signalling_models/KS_network.tab', sep='\t', header=0, dtype=str).dropna()
ks_network['RID'] = [x[0] for x in ks_network['SID']]
ks_network['site'] = ks_network['res'] + ks_network['pos']
print '[INFO] network: ', ks_network.shape

# Filter predicted interactions
ks_network = ks_network.loc[[not x.startswith('n') for x in ks_network['SID']]]
print '[INFO] network, remove predicted interactions: ', ks_network.shape

# Remove self phosphorylations
ks_network = ks_network[[r['S.AC'] != r['K.AC'] for i, r in ks_network.iterrows()]]
print '[INFO] network, remove self phosphorylations: ', ks_network.shape

# Calculate kinase target sites
kinases = set(ks_network.loc[ks_network['K.AC'] != '', 'K.AC'])
kinases_targets = {k: set(map('_'.join, ks_network.loc[ks_network['K.AC'] == k, ['S.AC', 'site']].values)).intersection(pp_fc.index) for k in kinases}
kinases_targets = {k: v for k, v in kinases_targets.iteritems() if len(v) > 0}
print '[INFO] Kinases targets: ', len(kinases_targets)

# Calculate kinase activity score
kinases_es = {k: gsea(pp_fc.to_dict(), kinases_targets[k], 10000) for k in kinases_targets}
kinases_es = {k: -np.log10(kinases_es[k][1]) if kinases_es[k][0] < 0 else np.log10(kinases_es[k][1]) for k in kinases_es}
print '[INFO] Kinase enrichment: ', len(kinases_es)

# ---- Perform GSEA on the kinase activities
# Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()
kegg.organism = 'hsa'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Convert KEGG pathways Gene Name to UniProt
kegg_to_uniprot = kegg.conv('uniprot', 'hsa')
kegg_pathways = {p: {kegg_to_uniprot[x].split(':')[1] for i in kegg_pathways[p]['entries'] if i['type'] == 'gene' for x in i['name'].split(' ') if x in kegg_to_uniprot} for p in kegg_pathways}
kegg_pathways = {p: kegg_pathways[p] for p in kegg_pathways if len(kegg_pathways[p].intersection(kinases_es)) > 2}
kegg_pathways_names = {p: re.findall('NAME\s*(.*) - Homo sapiens\n?', kegg.get(p))[0] for p in kegg_pathways}

# ---- Pathway enrichment analysis
pathways_es = {p: gsea(kinases_es, kegg_pathways[p], 10000) for p in kegg_pathways}
pathways_es = DataFrame(pathways_es, index=['es', 'pvalue']).T
pathways_es['adj.pvalue'] = multipletests(pathways_es['pvalue'], method='fdr_bh')[1]
pathways_es['name'] = [kegg_pathways_names[i] for i in pathways_es.index]
pathways_es['intersection'] = [len(kegg_pathways[i].intersection(kinases_es.keys())) for i in pathways_es.index]
pathways_es = pathways_es.sort('adj.pvalue')
print '[INFO] Pathways enrichment: ', len(pathways_es)
