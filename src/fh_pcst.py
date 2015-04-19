__author__ = 'emanuel'

import time
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
import pydot
from pandas import read_csv, Series
from pcst.network import Network, Type, Edge
from pcst.formulation import PCST
from utils.map_peptide_sequence import match_sequence, read_fasta, read_uniprot_genename, read_uniprot_accname


def hill_function(matrix, hill_coef=2):
    return 1 / ((np.median(matrix, axis=0) / matrix) ** hill_coef + 1)


def map_uniprot(uniprots):
    return {} if len(uniprots) == 0 else {uniprot2gene[u] for u in uniprots if u in uniprot2gene}

# Configure
data_wd = '/Users/emanuel/Projects/projects/pymist/reports/fh_cells/'
p_value_thres = 0.05
reac_thres = 1
inverse_penalty = np.Inf
uniprot_fasta = read_fasta()

# Import data-sets
pp_human = read_csv(data_wd + 'b1368p100_phospho_human_limma.tsv', sep='\t')
tp_human = read_csv(data_wd + 'b1368p100_protein_human_limma.tsv', sep='\t')

# Remove protein ambigous phosphosites
pp_human = pp_human[[len(x.split('; ')) <= 2 for x in pp_human['uniprot']]]

# Uniprot to gene name map
uniprots_sequences = read_fasta()
uniprot2gene = read_uniprot_genename()
accname2uniprot = dict((value, key) for key, value in read_uniprot_accname().iteritems())
genename2uniprot = dict((value, key) for key, value in read_uniprot_genename().iteritems())

# Annotate peptides with uniprot IDs
pp_human['uniprot'] = [x.split('; ')[0] for x in pp_human['uniprot']]
tp_human['uniprot'] = [accname2uniprot[x] for x in tp_human['acc_no']]

protein_measured = set(pp_human[pp_human['adj.P.Val'] < 0.05]['uniprot'])
protein_measured_gene = [uniprot2gene[p] for p in protein_measured if p in uniprot2gene]

# Metabolic genes
recon1_genes = read_csv('/Users/emanuel/Projects/resources/metabolic_models/recon1_genes.txt', sep='\t', names=['gene_symbol'])
recon1_genes = set(recon1_genes.gene_symbol)

# Import network
ks_network = read_csv('/Users/emanuel/Projects/resources/signalling_models/KS_network.tab', sep='\t', header=0, dtype=str)
ks_network = ks_network.dropna()
ks_network = ks_network.loc[[not x.startswith('n') for x in ks_network['SID']]]
ks_network['RID'] = [x[0] for x in ks_network['SID']]

# Initialise network model
network = Network()

# Build network
for index, row in ks_network.iterrows():
    kid = str(row['K.AC'])
    sid = str(row['S.AC'])
    res = (str(row['res']), int(row['pos']))
    database = str(row['RID'])

    # Discard self-phosphorylations
    if kid != sid:
        # Create kinase node
        if not network.has_node(kid):
            if kid in protein_measured:
                kid_value = np.max(np.abs(pp_human[pp_human['uniprot'] == kid]['logFC']))
                network.create_node(kid, kid_value, Type.TERMINAL)

            else:
                network.create_node(kid)

        # Create kubstrate node
        if not network.has_node(sid):
            if sid in protein_measured:
                sid_value = np.max(np.abs(pp_human[pp_human['uniprot'] == sid]['logFC']))
                network.create_node(sid, sid_value, Type.TERMINAL)

            else:
                network.create_node(sid)

        # Add interaction or increment value
        edge = network.get_edge(kid + '_' + sid)

        if (edge is not None) and (not database in edge.databases):
            # edge.value += 1
            edge.databases.append(database)

        else:
            network.create_edge(kid, sid, 0.5, databases=[database])

print 'KS network'
print network

# Normalise nodes values
nodes_penalties = Series({n.id: n.value for n in network.nodes.values() if n.value != 0.0})
nodes_penalties = hill_function(nodes_penalties)
for i in nodes_penalties.index:
    network.get_node(i).value = nodes_penalties[i]

# Normalise nodes values
edges_penalties = Series({e.id: e.value for e in network.edges.values() if e.value != 10})
for e in network.edges.values():
    e.value -= network.nodes[e.target].value

# Convert to directed graph
network.convert_to_directed(inverse_penalty)

print 'Inverse edges added'
print network

start = time.time()

# Formulate PCST
pcst = PCST(network)
pcst.initialise()

# Run PCST
beta = 1
pcst.objective_function(beta)
solution = pcst.optimise()

print 'LP solved! Elapsed time: ', str(int(time.time() - start)), 's'

# Percentage of network nodes that are measured
solution_nodes = {e.source for e in solution.edges.values()}.union({e.target for e in solution.edges.values()})
solution_nodes_measured = [n for n in solution_nodes if network.get_node(n).type == Type.TERMINAL]
print 'Network nodes measured: %.2f%%' % (1.0 * len(solution_nodes_measured) / len(solution_nodes) * 100)

# Percentage of measured nodes that are in the network
protein_measured_solution = [n for n in protein_measured if n in solution_nodes]
print 'Measured nodes in solution: %.2f%%' % (1.0 * len(protein_measured_solution) / len(protein_measured) * 100)

# Number of interactions
print 'Number of interactions: %d' % len(solution.edges)

# Metabolic enzymes
metabolic_targets = [uniprot2gene[n] for n in solution_nodes if n in uniprot2gene and uniprot2gene[n] in recon1_genes]

# Check back to experimental data
gene_name = 'ABL1'
gene_uniprot = [k for k, v in uniprot2gene.items() if v == gene_name][0]
pp_human.loc[[gene_uniprot in u for u in pp_human.uniprot]]

# Check interactions
gene_name = 'ABL1'
[uniprot2gene[e.target] for e in network.out_edges([k for k, v in uniprot2gene.items() if v == gene_name][0]) if e.target in uniprot2gene]

# # Edge sensitivity analysis
# edge_sensitivity_scores = pcst.edge_sensitivity_analysis(solution)

# Save network as pdf
solution.to_pdf(data_wd + 'pp_human_dp_network_' + 'beta=' + str(beta) + '.pdf', uniprot2gene, metabolic_targets, protein_measured_gene)


# def plot_interactions(pairs):
pairs = [('ABL1', 'PDK1'), ('PDK1', 'PDHA1')]

graph = pydot.Dot('graphname', graph_type='digraph', compound='true', rankdir='LR')

for kinase, substrate in pairs:
    k_uniprot, s_uniprot = genename2uniprot[kinase], genename2uniprot[substrate]
    interactions = ks_network[np.bitwise_and(ks_network['S.AC'] == s_uniprot, ks_network['K.AC'] == k_uniprot)]

    k_cluster = pydot.Cluster(k_uniprot, label=kinase, style='rounded', color='')
    graph.add_subgraph(k_cluster)

    dummy = pydot.Node(kinase + '_invis', style='invis', width='0.011', height='0.021', fixedsize='true')
    k_cluster.add_node(dummy)

    s_cluster = pydot.Cluster(s_uniprot, label=substrate, style='rounded')
    graph.add_subgraph(s_cluster)

    for site in set(interactions['res'].values + interactions['pos'].values):
        node = pydot.Node(site, style='filled')
        s_cluster.add_node(node)

        site_logfc = pp_human[[x.startswith(substrate + '(' + site + ')') for x in pp_human['site']]]

        if len(site_logfc) > 0 and substrate in metabolic_targets:
            node.set_fillcolor('#8EC127')
            node.set_color('#FFCC00')

        elif len(site_logfc) > 0 and substrate in protein_measured_gene:
            node.set_fillcolor('#3498db')
            node.set_color('#FFCC00')

        else:
            node.set_fillcolor('gray')

        graph.add_edge(pydot.Edge(dummy, node, ltail=k_cluster.get_name()))

graph.write_pdf('/Users/emanuel/Downloads/test.pdf')
