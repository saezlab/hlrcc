import re
import pydot
import igraph
import numpy as np
import seaborn as sns
import subprocess
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
# c_info = ['KIN_ACC_ID', 'SUB_ACC_ID', 'SUB_MOD_RSD', 'KIN_ORGANISM', 'SUB_ORGANISM']
#
# network = read_csv(wd + '/files/kinase_substrate_phosphositeplus.txt', sep='\t')[c_info]
# network['SITE'] = network['SUB_ACC_ID'] + '_' + network['SUB_MOD_RSD']
# print '[INFO] Kinase/substrate network: ', len(network)
#
# network = network[np.bitwise_and(network['KIN_ORGANISM'] == 'human', network['SUB_ORGANISM'] == 'human')]
# print '[INFO] Kinase/substrate network, human: ', len(network)
#
# network = network[network['KIN_ACC_ID'] != network['SUB_ACC_ID']]
# print '[INFO] Kinase/substrate network, human, no self-phosphorylations: ', len(network)
#
# vertices = list(set(network['KIN_ACC_ID']).union(network['SUB_ACC_ID']).union(network['SITE']))

network = read_csv('/Users/emanuel/Projects/resources/signalling_models/KS_network.tab', sep='\t', header=0, dtype=str).dropna()
network = network[[not x.startswith('n') for x in network['SID']]]
network = network[[s != k for s, k in network[['S.AC', 'K.AC']].values]]
network['SITE'] = network['S.AC'] + '_' + network['res'] + network['pos']
print '[INFO] network: ', network.shape

vertices = list(set(network['K.AC']).union(network['S.AC']).union(network['SITE']))

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
    # source, site, substrate = network.ix[i, 'KIN_ACC_ID'], network.ix[i, 'SITE'], network.ix[i, 'SUB_ACC_ID']
    source, site, substrate = network.ix[i, 'K.AC'], network.ix[i, 'SITE'], network.ix[i, 'S.AC']

    edges.append((source, site))
    edges.append((site, substrate))

network_i.add_edges(edges)
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print '[INFO] Sub-network largest component simplified: ', network_i.summary()

# # Get largest component
# network_i = network_i.components().giant()
# print '[INFO] Sub-network largest component: ', network_i.summary()

# Export network to Heinz
with open('%s/files/heinz_nodes.txt' % wd, 'w') as f:
    f.write('#node\tscore1\n')
    f.write('\n'.join(['%s\t%.4f' % (v['name'], v['weight']) for v in network_i.vs]))

with open('%s/files/heinz_edges.txt' % wd, 'w') as f:
    f.write('#nodeA\tnodeB\tscore1\n')
    f.write('\n'.join(['%s\t%s\t0.0' % tuple(network_i.vs[[e.source, e.target]]['name']) for e in network_i.es]))

print '[INFO] Network exported to Heinz'

# Run Heinz
subprocess.call('heinz -e %s/files/heinz_edges.txt -n %s/files/heinz_nodes.txt -o %s/files/heinz_solution.txt -t 3000' % (wd, wd, wd), shell=True)
print '[INFO] Heinz execution finished'

# Read Heinz solution
heinz_network = read_csv('%s/files/heinz_solution.txt' % wd, sep='\t')[:-1].dropna()
heinz_network = network_i.subgraph(heinz_network['#label'])
print '[INFO] Heinz active module: ', heinz_network.summary()

# Draw network
graph = pydot.Dot(graph_type='graph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3')
graph.set_edge_defaults(color='gray', arrowhead='vee')

for edge in heinz_network.es:
    source_id, target_id = heinz_network.vs[[edge.source, edge.target]]['name']

    source = pydot.Node(source_id, style='filled', shape='box', penwidth='0')
    target = pydot.Node(target_id, style='filled', shape='box', penwidth='0')

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('%s/files/heinz_solution.pdf' % wd)
print '[INFO] Network PDF saved!\n'

# [heinz_network.vs[x]['name'] for x in heinz_network.get_all_shortest_paths('P00519', 'P08559')]
# [network_i.vs[x]['name'] for x in network_i.get_all_shortest_paths('P00519', 'P08559')]

# ABL2: P42684, LASP1: Q14847
# P00519, P42684_Y261
# P00519, Q14847_Y171

print heinz_network.get_eid('P00519', 'P42684_Y261')
print heinz_network.get_eid('P00519', 'Q14847_Y171')