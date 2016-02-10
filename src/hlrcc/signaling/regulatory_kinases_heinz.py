import re
import pydot
import igraph
import pickle
import numpy as np
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from bioservices import UniProt
from pymist.reader.sbml_reader import read_sbml_model
from pymist.utils.omnipath_phospho import get_targets
from statsmodels.distributions import ECDF
from pandas import DataFrame, Series, read_csv


# -- Import data-sets
# Phospho fold-change
phospho_fc = read_csv('%s/data/uok262_phosphoproteomics_logfc.txt' % wd, sep='\t')
phospho_fc['psite'] = ['_'.join(i.split('_')[:2]) for i in phospho_fc.index]
phospho_fc = phospho_fc.groupby('psite')['logFC'].median()

# Import kinase activity
k_activity = read_csv('%s/data/uok262_kinases_activity_lm.txt' % wd, sep='\t', index_col=0, names=['kinase', 'activity'])

# Kinases targets
sources = ['HPRD', 'PhosphoSite', 'Signor', 'phosphoELM']
k_targets = get_targets(sources, remove_self=False)
print '[INFO] Kinases targets imported: ', k_targets.shape


# -- Build Kinase/Substrate network
network = k_targets.replace(0, np.nan)
network = [(k, v) for k in network for v in network[k].dropna().index]

vertices = {v for i in network for v in i}
vertices = list(vertices.union({v.split('_')[0] for v in vertices if len(v.split('_')) > 1}))


# -- Signalling network analysis
# Scale kinases activities and sites fold-changes
k_activity_abs, phospho_fc_abs = k_activity.abs(), phospho_fc.abs()

# Scale kinases enrichment
kinases_ecdf = ECDF(k_activity_abs['activity'])
kinases_weights = {k: kinases_ecdf(k_activity_abs.ix[k])[0] for k in k_activity_abs.index if k in vertices}

# Scale sites FC
phospho_ecdf = ECDF(phospho_fc_abs.values)
phospho_weights = {k: phospho_ecdf(phospho_fc_abs.ix[k]) for k in phospho_fc_abs.index if k in vertices}

# Vertices weights
vertices_weights = kinases_weights.copy()
vertices_weights.update(phospho_weights)


# -- Create network
network_i = igraph.Graph(directed=False)

# Add nodes
network_i.add_vertices(vertices)
network_i.vs['weight'] = [vertices_weights[v] if v in vertices_weights else -1.0 for v in vertices]
print '[INFO] Kinase/Substrate network: ', network_i.summary(), len(vertices_weights)

# Add edges
edges = []
for e in network:
    kinase, substrate = e

    edges.append((kinase, substrate))
    edges.append((substrate, substrate.split('_')[0]))

network_i.add_edges(edges)
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print '[INFO] Sub-network largest component simplified: ', network_i.summary()

# Get largest component
network_i = network_i.components().giant()
print '[INFO] Sub-network largest component: ', network_i.summary()

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

    for node in [source, target]:
        if node.get_name() in phospho_fc.index:
            node.set_fillcolor('#3498db')

        elif node.get_name() in k_activity_abs.index:
            node.set_fillcolor('#BB3011')

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('%s/reports/heinz_solution.pdf' % wd)
print '[INFO] Network PDF saved!\n'

# Sub-network
sub_network = network_i.subgraph(network_i.neighborhood(vertices=['P08559_S232'], order=3)[0])
print '[INFO] Subnetwork: ', sub_network.summary()

# [heinz_network.vs[x]['name'] for x in heinz_network.get_all_shortest_paths('P00519', 'P08559_S232')]
# [network_i.vs[x]['name'] for x in network_i.get_all_shortest_paths('P00519', 'P08559_S232')]

# ABL2: P42684, LASP1: Q14847
# P00519, P42684_Y261
# P00519, Q14847_Y171

print heinz_network.get_eid('P00519', 'P42684_Y261')
print heinz_network.get_eid('P00519', 'Q14847_Y171')
