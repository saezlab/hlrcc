import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pandas import DataFrame, Series, read_csv
from bioservices import KEGG, KEGGParser

sns.set_style('white')

# ---- Import data-sets
pp = read_csv('%s/files/b1368p100_phospho_human_limma.tsv' % wd, sep='\t')
tp = read_csv('%s/files/b1368p100_protein_human_limma.tsv' % wd, sep='\t')

# ---- Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()

kegg.organism = 'hsa'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Convert KEGG pathways Gene Name to UniProt
kegg_to_uniprot = kegg.conv('uniprot', 'hsa')
kegg_pathways_proteins = {p: {kegg_to_uniprot[x].split(':')[1] for i in kegg_pathways[p]['entries'] if i['type'] == 'gene' for x in i['name'].split(' ') if x in kegg_to_uniprot} for p in kegg_pathways}
print '[INFO] KEGG genes converted to UniProt: ', len(kegg_pathways_proteins)

