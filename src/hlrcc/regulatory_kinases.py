import re
import pydot
import igraph
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from hlrcc.figures.volcano_plots import get_fasta
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

# ---- Import metabolic enzymes
enzymes = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml').get_genes()

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

# ---- Import uniprot Bioservice
uniprot_bioservice = UniProt(cache=True)


def uniprot2genename(uniprotid):
    try:
        res = re.findall('.* GN=(.*?) ', uniprot_bioservice.get_fasta(uniprotid))

    except TypeError:
        res = []

    return res[0] if len(res) == 1 else ''

print '[INFO] UniProt Bioservice established'

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

# ---- Create network
network_i = igraph.Graph(directed=True)

# Add nodes
network_i.add_vertices(vertices)
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Add edges
edges, edges_names, edges_weights = [], [], []
for i in network.index:
    source, site, substrate = network.ix[i, 'KIN_ACC_ID'], network.ix[i, 'SITE'], network.ix[i, 'SUB_ACC_ID']

    edges.append((source, site))
    edges_names.append('%s -> %s' % (source, site))
    edges_weights.append(1.0 - kinases_weights[source] if source in kinases_weights else 1.0)

    edges.append((site, substrate))
    edges_names.append('%s -> %s' % (site, substrate))
    edges_weights.append(1.0 - phospho_weights[site] if site in phospho_weights else 1.0)

network_i.add_edges(edges)
network_i.es['name'] = edges_names
network_i.es['weight'] = edges_weights
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# Remove duplicated edges and self-loops
network_i = network_i.simplify(True, True, 'first')
print '[INFO] Kinase/Substrate network: ', network_i.summary()

# ---- Sub-set network to differentially phosphorylated sites
sub_network = network_i.subgraph({x for i in phospho_fc[phospho_fc.abs() > .5].index if i in vertices for x in network_i.neighborhood(i, order=5, mode='IN')})
print '[INFO] Sub-network created: ', sub_network.summary()

sub_network_nweighted = sub_network.copy()
print '[INFO] Unweighted sub-network created: ', sub_network_nweighted.summary()

sub_network_weighted = sub_network.spanning_tree('weight')
print '[INFO] Weighted sub-network created: ', sub_network_weighted.summary()


# ---- Plot consensus network
def draw_network(network2draw, path):
    graph = pydot.Dot(graph_type='digraph', rankdir='LR')

    graph.set_node_defaults(fontcolor='white', penwidth='3')
    graph.set_edge_defaults(color='gray', arrowhead='vee')

    freq_ecdf = ECDF(network2draw.es.get_attribute_values('weight'))

    for edge_index in network2draw.es.indices:
        edge = network2draw.es[edge_index]

        source_id, target_id = network2draw.vs[edge.source].attributes()['name'], network2draw.vs[edge.target].attributes()['name']

        source = pydot.Node(source_id, style='filled', shape='box', penwidth='0')
        target = pydot.Node(target_id, style='filled', shape='box', penwidth='0')

        for node in [source, target]:
            node_name = node.get_name()
            protein_name = uniprot2genename(node_name.split('_')[0])

            # Set node colour
            if protein_name in enzymes:
                node.set_fillcolor('#8EC127')

            elif node_name in phospho_fc.index:
                node.set_fillcolor('#3498db')

            elif node_name in kinases_es:
                node.set_fillcolor('#BB3011')

            # Set node name
            if len(node_name.split('_')) == 2:
                node_name, res = node_name.split('_')
                node.set_name((protein_name if protein_name != '' else node_name) + '_' + res)

            else:
                node.set_name(protein_name if protein_name != '' else node_name)

            graph.add_node(node)

        # Set edge width
        edge_width = str((1 - freq_ecdf(network2draw.es[edge_index].attributes()['weight'])) * 5 + 1)

        edge = pydot.Edge(source, target, penwidth=edge_width)
        graph.add_edge(edge)

    graph.write_pdf(path)
    print '[INFO] Network PDF saved!\n'

draw_network(sub_network_weighted, wd + '/reports/signalling_regulatory_network.pdf')