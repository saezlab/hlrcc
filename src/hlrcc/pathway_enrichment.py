import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from statsmodels.stats.multitest import multipletests
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv
from bioservices import KEGG, KEGGParser, UniProt

# ---- Import Kinase enrichment
kinases_es = read_csv(wd + '/files/kinase_activity.tab', sep='\t')

# ---- Perform GSEA on the kinase activities
# Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()
kegg.organism = 'hsa'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Set-up UniProt bioservice
uniprot = UniProt(cache=True)

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

p = 'path:hsa05211'
gsea(kinases_es, kegg_pathways[p], permutations=1000, plot_name='%s/reports/kegg_pathway_kinase_activity_enrichment.pdf' % wd, plot_title=kegg_pathways_names[p])

for protein in kegg_pathways[p]:
    print protein
    print re.findall('.*\|(.*)', uniprot.get_fasta(protein).split(' ')[0])[0]