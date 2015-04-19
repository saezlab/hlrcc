import re
from pandas.stats.misc import zscore
import pydot
import numpy as np
import igraph as igraph
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
from pandas import DataFrame, Series, read_csv, pivot_table
from enrichment.gsea import gsea
from utils.map_peptide_sequence import read_fasta, read_uniprot_genename


wd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/'
data_wd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/files/'

# Import Uniprot
uniprot2gene = read_uniprot_genename()

# Metabolic enzymes
metabolic_targets = set(read_csv('/Users/emanuel/Projects/resources/metabolic_models/recon1_genes.txt', sep='\t', names=['gene_symbol'])['gene_symbol'])
metabolic_targets = [uniprot for uniprot, (gene_name, acc_number) in uniprot2gene.items() for gene in metabolic_targets if gene_name == gene]

# Import phospho
pp_human = read_csv(data_wd + 'b1368p100_phospho_human_limma.tsv', sep='\t')
print '[INFO] pp: ', pp_human.shape

# Remove ambigous mapping
pp_human = pp_human[[len(s.split(';')) == 2 for s in pp_human['site']]]
print '[INFO] pp, ambigous mapping: ', pp_human.shape

# Remove multiple-phosphorylation
pp_human = pp_human[[len(s.split('+')) == 1 for s in pp_human['site']]]
print '[INFO] pp, multiple-phosphorylation: ', pp_human.shape

# Import total protein
tp_human = read_csv(data_wd + 'b1368p100_protein_human_limma.tsv', sep='\t')
print '[INFO] tp: ', tp_human.shape

# Remove ambigous mapping
tp_human = tp_human[[len(str(s).split(';')) == 2 for s in tp_human['uniprot']]]
tp_human['uniprot'] = [s.split(';')[0].strip() for s in tp_human['uniprot']]
print '[INFO] tp, ambigous mapping: ', tp_human.shape

# Average protein concentration
tp_fc = tp_human.groupby('uniprot')['logFC'].median()[1:]
print '[INFO] tp, average peptides: ', tp_fc.shape

# Parse sites
pp_human['site'] = [re.findall('\(.*\)', s)[0][1:].replace(')', '') for s in pp_human['site']]
pp_human['uniprot'] = [s.split(';')[0].strip() for s in pp_human['uniprot']]
pp_human['site'] = pp_human['uniprot'] + '_' + pp_human['site']
pp_fc = pp_human.groupby('site')['logFC'].median()
print '[INFO] pp, average phosphosites: ', pp_fc.shape

pp_fc = Series(dict(zip(pp_fc.index, pp_fc.ix[np.random.randint(0, len(pp_fc), len(pp_fc))].values)))

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
kinases_targets = {k: v for k, v in kinases_targets.iteritems() if len(v) > 2}

# Parse columns
ks_network['site'] = ks_network['S.AC'] + '_' + ks_network['site']
ks_network['kinase'] = ks_network['K.AC']
ks_network['substrate'] = ks_network['S.AC']

# Calculate kinase change
kinases_es = {k: -np.log10(gsea(pp_fc, targets, True, 1000)[1]) for k, targets in kinases_targets.items()}
print '[INFO] Kinase enrichment: ', len(kinases_es)

# Scale fold-changes
sites_fc_ecdf = ECDF(pp_fc.values)
sites_fc = {k: sites_fc_ecdf(v) for k, v in pp_fc.to_dict().items()}

# Scale kinase enrichments
kinases_es_ecdf = ECDF(kinases_es.values())
kinases_es = {k: kinases_es_ecdf(v) for k, v in kinases_es.items()}

# Create network
vertices = set(ks_network['kinase']).union(set(ks_network['substrate'])).union(set(ks_network['site']))
vertices = dict(zip(*(vertices, range(len(vertices)))))
vertices_inv = {v: k for k, v in vertices.items()}

edges, edges_weights = [], []
for i, row in ks_network.iterrows():
    edges.append((vertices[row['kinase']], vertices[row['site']]))
    edges_weights.append(1.0 - kinases_es[row['kinase']] if row['kinase'] in kinases_es else 1.0)

    edges.append((vertices[row['site']], vertices[row['substrate']]))
    edges_weights.append(1.0 - sites_fc[row['site']] if row['site'] in sites_fc else 1.0)

network = igraph.Graph(edges, directed=True, edge_attrs={'weight': edges_weights})
network.simplify(True, False, 'first')
network.vs.set_attribute_values('name', [vertices_inv[k.index] for k in network.vs])
print '[INFO] Network created: ', network.summary()

kinases = set(ks_network.loc[ks_network['K.AC'] != '', 'K.AC'])
sites = set(ks_network['site']).intersection(pp_fc.index)

shortest_paths = [(kinase, site, network.get_all_shortest_paths(vertices[kinase], vertices[site], 'weight')) for kinase in kinases for site in sites]
print '[INFO] Shortest paths calculated: ', len(shortest_paths)

shortest_paths_all = [(kinase, site, network.get_all_shortest_paths(vertices[kinase], vertices[site])) for kinase in kinases for site in sites]
print '[INFO] Shortest paths calculated: ', len(shortest_paths_all)

shortest_paths_len = [(kinase, site, network.shortest_paths(vertices[kinase], vertices[site], 'weight')) for kinase in kinases for site in sites]
print '[INFO] Lenght of shortest paths calculated: ', len(shortest_paths_len)

# Consensus network analysis
shortest_paths_edges = [p for k, s, paths in shortest_paths for path in paths for p in zip(path, path[1:])]
shortest_paths_edges = [network.get_eid(s, t) for s, t in shortest_paths_edges]
shortest_paths_edges_freq = np.unique(shortest_paths_edges, return_counts=True)
shortest_paths_edges_freq = dict(zip(shortest_paths_edges_freq[0], np.log2(shortest_paths_edges_freq[1])))

# Unweighted shortest paths
shortest_paths_edges_all = [p for k, s, paths in shortest_paths_all for path in paths for p in zip(path, path[1:])]
shortest_paths_edges_all = [network.get_eid(s, t) for s, t in shortest_paths_edges_all]
shortest_paths_edges_all_freq = np.unique(shortest_paths_edges_all, return_counts=True)
shortest_paths_edges_all_freq = dict(zip(shortest_paths_edges_all_freq[0], np.log2(shortest_paths_edges_all_freq[1])))

# Save edges frequency as attribute
network.es.set_attribute_values('freq', [(1.0 * shortest_paths_edges_freq[index]) if index in shortest_paths_edges_freq else np.NaN for index in network.es.indices])

sns.set(style='white')
sns.distplot(shortest_paths_edges_freq.values(), color='#626266')
sns.despine()
plt.ylabel('log2 (frequency of edges counts)')
plt.xlabel('edges counts')
plt.savefig(wd + 'plots/shortest_paths_edges_freq.pdf', bbox_inches='tight')
plt.close('all')

sns.set(style='white')
plt.plot(range(16), [np.sum([v > cutoff for v in shortest_paths_edges_freq.values()]) for cutoff in range(16)], '-o', c='#626266', lw=2, ls='--')
sns.despine()
plt.ylabel('number of edges')
plt.xlabel('cut-off (log2 edges counts)')
plt.savefig(wd + 'plots/shortest_paths_edges_cutoff.pdf', bbox_inches='tight')
plt.close('all')

sns.set(style='white')
x = range(16)
y = [network.subgraph_edges([k for k, v in shortest_paths_edges_freq.items() if v > cutoff]).density() for cutoff in x]
plt.plot(x, y, '-o', c='#626266', lw=2, ls='--')
sns.despine()
plt.ylabel('network density')
plt.xlabel('cut-off (log2 edges counts)')
plt.savefig(wd + 'plots/shortest_paths_edges_cutoff_density.pdf', bbox_inches='tight')
plt.close('all')

sns.set(style='white')
x = range(16)
y = [len(network.subgraph_edges([k for k, v in shortest_paths_edges_freq.items() if v > cutoff]).clusters(mode='week')) for cutoff in x]
plt.plot(x, y, '-o', c='#626266', lw=2, ls='--')
sns.despine()
plt.ylabel('network clusters (week)')
plt.xlabel('cut-off (log2 edges counts)')
plt.savefig(wd + 'plots/shortest_paths_edges_cutoff_clusters.pdf', bbox_inches='tight')
plt.close('all')

# Shortest paths length analysis
shortest_paths_len_df = DataFrame([(k, s, l[0][0]) for k, s, l in shortest_paths_len if l[0][0] != np.Inf], columns=['kinase', 'site', 'length'])
shortest_paths_len_df = pivot_table(shortest_paths_len_df, values='length', index='kinase', columns='site')

sns.clustermap(shortest_paths_len_df.T.replace(np.NaN, 0), figsize=(40, 35), linewidths=0)
plt.savefig(wd + 'plots/kinase_site_length.png', bbox_inches='tight')
plt.close('all')

# Calculate consensus network
cutoff = 10
sub_network = network.subgraph_edges([k for k, v in shortest_paths_edges_freq.items() if v > cutoff])
print '[INFO] Consensus cutoff: ', cutoff
print '[INFO] Consensus network: ', sub_network.summary()

# Plot
graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3')
graph.set_edge_defaults(color='gray', arrowhead='vee')

freq_ecdf = ECDF(sub_network.es.get_attribute_values('freq'))

for edge_index in sub_network.es.indices:
    edge = sub_network.es[edge_index]

    source_id, target_id = sub_network.vs[edge.source].attributes()['name'], sub_network.vs[edge.target].attributes()['name']

    source = pydot.Node(source_id, style='filled', shape='box', penwidth='0')
    target = pydot.Node(target_id, style='filled')

    for node in [source, target]:
        node_name = node.get_name()

        # Set node colour
        if node_name.split('_')[0] in metabolic_targets:
            node.set_fillcolor('#8EC127')

        elif node_name in sites_fc:
            node.set_fillcolor('#3498db')

        elif node_name in kinases_es:
            node.set_fillcolor('#BB3011')

        # Set node name
        if len(node_name.split('_')) == 2:
            node_name, res = node_name.split('_')
            node.set_name((uniprot2gene[node_name][0] if node_name in uniprot2gene else node_name) + '_' + res)

        else:
            node.set_name(uniprot2gene[node_name][0] if node_name in uniprot2gene else node_name)

        graph.add_node(node)

    # Set edge width
    edge_width = str(freq_ecdf(sub_network.es[edge_index].attributes()['freq']) * 5 + 1)

    edge = pydot.Edge(source, target, penwidth=edge_width)
    graph.add_edge(edge)

graph.write_pdf(wd + 'plots/consensus_network_' + str(cutoff) + '.pdf')
print '[INFO] Network PDF saved!'

# ABL1: P00519, PDK1: Q15118, PDHA1: P08559, PDPK1: O15530, ABL2: P42684, LASP1: Q14847
edge = network.get_eid(network.vs.find(name='Q15118'), network.vs.find(name='P08559_S232'))

[[network.vs[p].attributes()['name'] for p in path] for path in network.get_all_shortest_paths('P00519', 'P08559_S232')]
[[network.vs[p].attributes()['name'] for p in path] for path in network.get_all_shortest_paths('P00519', 'P08559_S232', 'weight')]

s = network.vs.find(name='P08559_S232').index
predecessors_edges_freq = [(p.attributes()['name'], int(np.exp2(sp_dict[network.get_eid(p, s)]))) for p in network.vs[s].predecessors() for sp_dict in [shortest_paths_edges_freq, shortest_paths_edges_all_freq]]
predecessors_edges_freq = DataFrame([(k, v) for k, v in predecessors_edges_freq if v > 1])