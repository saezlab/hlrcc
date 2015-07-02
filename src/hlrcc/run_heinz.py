import re
import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from bioservices import UniProt
from pymist.reader.sbml_reader import read_sbml_model
from statsmodels.distributions import ECDF
from pandas import DataFrame, Series, read_csv

sns.set_style('white')

# -- Import data-set and networks
# ---- Import phospho FC
phospho_fc = read_csv(wd + '/data/b1368p100_phospho_human_limma.tsv', sep='\t', index_col=0)['logFC']

# ---- Import kinase activity
kinases_es = read_csv(wd + '/files/kinase_activity.tab', sep='\t', index_col=0, header=None, names=['activity'])['activity']

# ---- Import Kinase/Substrate network
c_info = ['KIN_ACC_ID', 'SUB_ACC_ID', 'SUB_MOD_RSD', 'KIN_ORGANISM', 'SUB_ORGANISM']

network = read_csv(wd + '/files/kinase_substrate_phosphositeplus.txt', sep='\t')[c_info]
network['SITE'] = network['SUB_ACC_ID'] + '_' + network['SUB_MOD_RSD']
print '[INFO] Kinase/substrate network: ', len(network)

network = network[np.bitwise_and(network['KIN_ORGANISM'] == 'human', network['SUB_ORGANISM'] == 'human')]
print '[INFO] Kinase/substrate network, human: ', len(network)

network = network[network['KIN_ACC_ID'] != network['SUB_ACC_ID']]
print '[INFO] Kinase/substrate network, human, no self-phosphorylations: ', len(network)

vertices = list(set(network['KIN_ACC_ID']).union(network['SUB_ACC_ID']).union(network['SITE']))

# -- Signalling network analysis
# ---- Scale kinases activities and sites fold-changes
# Abs sites FC and kinases activities
kinases_es_abs, phospho_fc_abs = kinases_es.abs(), phospho_fc.abs()

# Scale kinases enrichment
kinases_ecdf = ECDF(kinases_es_abs.values)
kinases_weights = {k: kinases_ecdf(kinases_es_abs.ix[k]) for k in kinases_es_abs.index if k in vertices}

# Scale sites FC
phospho_ecdf = ECDF(phospho_fc_abs.values)
phospho_weights = {k: phospho_ecdf(phospho_fc_abs.ix[k]) for k in phospho_fc_abs.index if k in vertices}

# Vertices weights
vertices_weights = kinases_weights.copy()
vertices_weights.update(phospho_weights)

# ---- Create network
network_i = igraph.Graph(directed=False)

# Add nodes
network_i.add_vertices(vertices)
network_i.vs['weight'] = [vertices_weights[v] if v in vertices_weights else -1.0 for v in vertices]
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Add edges
edges = []
for i in network.index:
    source, site, substrate = network.ix[i, 'KIN_ACC_ID'], network.ix[i, 'SITE'], network.ix[i, 'SUB_ACC_ID']

    edges.append((source, site))
    edges.append((site, substrate))

network_i.add_edges(edges)
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print '[INFO] Sub-network largest component simplified: ', network_i.summary()

# # Sub-network with measured nodes
# subnetwork = network_i.subgraph(vertices_weights)
# print '[INFO] Sub-network: ', subnetwork.summary()

# Get largest component
subnetwork = network_i.components().giant()
print '[INFO] Sub-network largest component: ', subnetwork.summary()

# [subnetwork.vs[x]['name'] for x in subnetwork.get_all_shortest_paths('P00519', 'P08559')]

# Export network to Heinz
with open('%s/files/heinz_nodes.txt' % wd, 'w') as f:
    f.write('#node\tscore1\n')
    f.write('\n'.join(['%s\t%.4f' % (v['name'], v['weight']) for v in subnetwork.vs]))

with open('%s/files/heinz_edges.txt' % wd, 'w') as f:
    f.write('#nodeA\tnodeB\tscore1\n')
    f.write('\n'.join(['%s\t%s\t0.0' % (network_i.vs[e.source]['name'], network_i.vs[e.target]['name']) for e in subnetwork.es]))

print '[INFO] Network exported to Heinz'

# 'heinz -e ', edges_file, ' -n ', nodes_file, ' -o ', result_file, ' -t 1500'