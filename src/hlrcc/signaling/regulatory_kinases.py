#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pymist.utils.omnipath_phospho import get_targets
from statsmodels.distributions import ECDF
from pandas import DataFrame, read_csv


# -- Import Uniprot id mapping
umap = read_csv('./files/protein-coding_gene.txt', sep='\t').dropna(subset=['uniprot_ids'])
umap = umap.groupby('uniprot_ids')['symbol'].agg(lambda x: ';'.join([g for i in x for g in i.split('|')]))


# -- Import data-sets
# Phospho fold-change
phospho_fc = read_csv('./data/uok262_phosphoproteomics_logfc.txt', sep='\t')
phospho_fc['psite'] = ['_'.join(i.split('_')[:2]) for i in phospho_fc.index]
phospho_fc = phospho_fc.groupby('psite')['logFC'].median()

# Import kinase activity
k_activity = read_csv('./data/uok262_kinases_activity_gsea.txt', sep='\t', index_col=0, names=['kinase', 'activity'])

# Kinases targets
sources = ['HPRD', 'PhosphoSite', 'Signor', 'phosphoELM', 'DEPOD']
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
psites = ['P04406_S83']  # Q01581_S495, P04406_S83, P08559_S232

sub_network = network_i.subgraph(network_i.neighborhood(vertices=psites, order=3, mode='IN')[0])
print '[INFO] Subnetwork: ', sub_network.summary()

# Draw network
graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3', fillcolor='#CCCCCC')
graph.set_edge_defaults(color='#CCCCCC', arrowhead='vee')

for edge in sub_network.es:
    source_id, target_id = sub_network.vs[[edge.source, edge.target]]['name']

    source = pydot.Node(source_id, style='filled', shape='ellipse', penwidth='0')
    source.set_label(source_id.replace(source_id.split('_')[0], umap[source_id.split('_')[0]]))

    target = pydot.Node(target_id, style='filled', shape='ellipse', penwidth='0')
    target.set_label(target_id.replace(target_id.split('_')[0], umap[target_id.split('_')[0]]))

    for node in [source, target]:
        if node.get_name() in phospho_fc.index:
            node.set_fillcolor('#3498db')

        elif node.get_name() in k_activity_abs.index:
            node.set_fillcolor('#34495e')

    graph.add_node(source)
    graph.add_node(target)

    edge = pydot.Edge(source, target)
    graph.add_edge(edge)

graph.write_pdf('./reports/psite_upstream_pathway_%s.pdf' % '_'.join(psites))
print '[INFO] Network PDF saved!\n'

# Plot kinases activities
k_level1 = set(sub_network.subgraph(sub_network.neighborhood(vertices=psites, order=1, mode='IN')[0]).vs['name'])

k_level2 = set(sub_network.subgraph(sub_network.neighborhood(vertices=psites, order=3, mode='IN')[0]).vs['name']) - k_level1
k_level2 = {k for k in k_level2 if len(k.split('_')) == 1}

plot_df = [(k, 'level 1' if k in k_level1 else 'level 2', umap[k], k_activity.ix[k].values[0]) for k in sub_network.vs['name'] if len(k.split('_')) == 1 and k in k_activity.index]
plot_df = DataFrame(plot_df, columns=['uniprot', 'level', 'name', 'activity'])
plot_df = plot_df.sort('activity', ascending=False)
plot_df['name'] = ['PDK' if i.startswith('PDK') else i for i in plot_df['name']]
plot_df['name'] = ['PDP' if i.startswith('PDP') and not i.startswith('PDPK') else i for i in plot_df['name']]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, col='level', sharey=False, sharex=False, aspect=0.75, col_order=['level 1', 'level 2'])
g.map(sns.barplot, 'name', 'activity', color='#34495e', lw=0, ci=None)
g.map(plt.axhline, y=0, lw=.3, ls='-', alpha=0.7, color='gray')
g.set_titles('Kinase {col_name}')
sns.despine()
plt.savefig('./reports/psite_upstream_pathway_%s_k_activity.pdf' % '_'.join(psites), bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
