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
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Import Uniprot id mapping
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)

# -- Import metabolic model
m_genes = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml').get_genes()


# -- Import data-sets
# Phospho fold-change
phospho_fc = read_csv('%s/data/uok262_phosphoproteomics_logfc.txt' % wd, sep='\t')
phospho_fc['psite'] = ['_'.join(i.split('_')[:2]) for i in phospho_fc.index]
phospho_fc = phospho_fc.groupby('psite')['logFC'].median()

# Import kinase activity
k_activity = read_csv('%s/data/uok262_kinases_activity_gsea.txt' % wd, sep='\t', index_col=0, names=['kinase', 'activity'])

# Kinases targets
sources = ['HPRD', 'PhosphoSite', 'Signor', 'phosphoELM']
k_targets = get_targets(sources, remove_self=False)
print '[INFO] Kinases targets imported: ', k_targets.shape

# Enzymes with regulatory p-sites
phospho_fc_enzymes = phospho_fc.ix[[i for i in phospho_fc.index if i.split('_')[0] in human_uniprot and human_uniprot[i.split('_')[0]][0] in m_genes]]
phospho_fc_enzymes = phospho_fc_enzymes.ix[[i for i in phospho_fc_enzymes.index if i in k_targets.index]].sort(inplace=False)


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
network_i = igraph.Graph(directed=True)

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


# -- Analyse p-site up-stream signalling pathway
psites = ['P08559_S232']  # Q01581_S495

sub_network = network_i.subgraph(network_i.neighborhood(vertices=psites, order=3, mode='IN')[0])
print '[INFO] Subnetwork: ', sub_network.summary()

# Draw network
graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3')
graph.set_edge_defaults(color='gray', arrowhead='vee')

for edge in sub_network.es:
    source_id, target_id = sub_network.vs[[edge.source, edge.target]]['name']

    source = pydot.Node(source_id, style='filled', shape='ellipse', penwidth='0')
    target = pydot.Node(target_id, style='filled', shape='ellipse', penwidth='0')

    for node in [source, target]:
        if node.get_name() in phospho_fc.index:
            node.set_fillcolor('#3498db')

        elif node.get_name() in k_activity_abs.index:
            node.set_fillcolor('#BB3011')

        n_name = node.get_name()
        n_name = n_name.replace(n_name.split('_')[0], human_uniprot[n_name.split('_')[0]][0])

        node.set_name(n_name)

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('%s/reports/psite_upstream_pathway_%s.pdf' % (wd, '_'.join(psites)))
print '[INFO] Network PDF saved!\n'

# Plot kinases activities
k_level1 = set(sub_network.subgraph(sub_network.neighborhood(vertices=psites, order=1, mode='IN')[0]).vs['name'])

k_level2 = set(sub_network.subgraph(sub_network.neighborhood(vertices=psites, order=3, mode='IN')[0]).vs['name']) - k_level1
k_level2 = {k for k in k_level2 if len(k.split('_')) == 1}

plot_df = [(k, 'level 1' if k in k_level1 else 'level 2', human_uniprot[k][0], k_activity.ix[k].values[0]) for k in sub_network.vs['name'] if len(k.split('_')) == 1 and k in k_activity.index]
plot_df = DataFrame(plot_df, columns=['uniprot', 'level', 'name', 'activity'])
plot_df = plot_df.sort('activity', ascending=False)

sns.set(style='ticks')
sns.barplot('name', 'activity', 'level', plot_df, palette=sns.light_palette('#34495e', 3)[1:], lw=0)
plt.axhline(0, lw=.3, ls='-', alpha=0.7, color='gray')
sns.despine()
plt.savefig('%s/reports/psite_upstream_pathway_%s_k_activity.pdf' % (wd, '_'.join(psites)), bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
