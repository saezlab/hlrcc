import re
import pydot
import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from gurobipy.gurobipy import GRB, LinExpr, QuadExpr, Model
from statsmodels.distributions.empirical_distribution import ECDF
from pcst.network import Network, Type, Edge
from pandas import DataFrame, Series, read_csv
from enrichment.gsea import gsea
from utils.map_peptide_sequence import read_fasta, read_uniprot_genename


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

# Parse sites
pp_human['site'] = [re.findall('\(.*\)', s)[0][1:].replace(')', '') for s in pp_human['site']]
pp_human['uniprot'] = [s.split(';')[0].strip() for s in pp_human['uniprot']]
pp_human['site'] = pp_human['uniprot'] + '_' + pp_human['site']
pp_fc = pp_human.groupby('site')['logFC'].median()
print '[INFO] pp, average phosphosites: ', pp_fc.shape

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

# Correct by protein concentration
# pp_fc = Series({k: (v / tp_fc.ix[k.split('_')[0]]) if k.split('_')[0] in tp_fc else v for k, v in pp_fc.to_dict().items()})

# Calculate kinase target sites
kinases = set(ks_network.loc[ks_network['K.AC'] != '', 'K.AC'])
kinases_targets = {k: set(map('_'.join, ks_network.loc[ks_network['K.AC'] == k, ['S.AC', 'site']].values)).intersection(pp_fc.index) for k in kinases}
kinases_targets = {k: v for k, v in kinases_targets.iteritems() if len(v) > 2}

# Calculate kinase change
kinases_es = {k: -np.log10(gsea(pp_fc, targets, True, 1000)[1]) for k, targets in kinases_targets.items()}
print '[INFO] Kinase enrichment: ', len(kinases_es)

ks_network['site'] = ks_network['S.AC'] + '_' + ks_network['site']
ks_network['kinase'] = ks_network['K.AC']
ks_network['substrate'] = ks_network['S.AC']

# Scale fold-changes
sites_fc_ecdf = ECDF(pp_fc.values)
sites_fc = {k: sites_fc_ecdf(v) for k, v in pp_fc.to_dict().items()}

# Scale kinase enrichments
kinases_es_ecdf = ECDF(kinases_es.values())
kinases_es = {k: kinases_es_ecdf(v) for k, v in kinases_es.items()}

# Create network
network = Network()
for i, row in ks_network.iterrows():
    kid, siteid, sid = row['kinase'], row['site'], row['substrate']

    kinase_es = kinases_es[kid] if kid in kinases_es else 0.0
    site_fc = sites_fc[siteid] if siteid in sites_fc else 0.0

    network.create_node(kid, kinase_es, 'protein')
    network.create_node(siteid, site_fc, 'site')
    network.create_node(sid, 0.0, 'protein')

    network.create_edge(kid, siteid, 1 - kinase_es)
    network.create_edge(siteid, sid, 1 - site_fc)

print '[INFO] Network created: \n' + str(network)

# Create problem
lp = Model('KS')

# Initialise nodes variables
for k, v in network.nodes.items():
    v.var = lp.addVar(vtype=GRB.BINARY, name=k)

for k, v in network.edges.items():
    v.var = lp.addVar(vtype=GRB.BINARY, name=k)

# Update variables
lp.update()
print '[INFO] Variables created!'

# Initialise constraints
for k, n in network.nodes.items():
    if n.type == 'site':
        lp.addConstr(LinExpr([(1, network.get_edge(edge).var) for edge in n.in_edges]), GRB.EQUAL, n.var)
        lp.addConstr(LinExpr([(1, network.get_edge(edge).var) for edge in n.out_edges]), GRB.EQUAL, n.var)

    elif n.type == 'protein':
        protein_targets, protein_sites = n.out_edges, n.in_edges

        if len(protein_sites) > 0:
            [lp.addConstr(network.get_edge(edge).var, GRB.LESS_EQUAL, n.var) for edge in protein_sites]
            lp.addConstr(LinExpr([(1, network.get_edge(edge).var) for edge in protein_sites]), GRB.GREATER_EQUAL, n.var)

        # if len(protein_targets) > 0:
        #     lhs = [(1, network.get_edge(edge).var) for edge in protein_targets]
        #     lp.addConstr(LinExpr(len(lhs), n.var), GRB.GREATER_EQUAL, LinExpr(lhs))

print '[INFO] Constraints created!'

# Initialise objective function
beta = 4
obj = LinExpr()

for edge in network.edges.values():
    obj.addTerms(edge.value, edge.var)

for node in network.nodes.values():
    if node.type == 'site':
        obj.addTerms(-beta * node.value, node.var)

lp.setObjective(obj, GRB.MINIMIZE)

# Optimise
lp.optimize()

solution = {var.VarName: var.x for var in lp.getVars() if var.x != 0.0}
solution_all = {var.VarName: var.x for var in lp.getVars()}
print '[INFO] Optimisation done: %.2f%% measured nodes' % (100.0 * len(set(solution).intersection(sites_fc)) / len(sites_fc))

graph = pydot.Dot(graph_type='digraph', rankdir='LR')

graph.set_node_defaults(fontcolor='white', penwidth='3')
graph.set_edge_defaults(color='gray', arrowhead='vee')

for e, v in solution.items():
    if len(e.split('->')) > 1:
        source_id, target_id = e.split('->')[0], e.split('->')[1]

        source = pydot.Node(source_id, style='filled', shape='box', penwidth='0')
        target = pydot.Node(target_id, style='filled')

        for node in [source, target]:
            if node.get_name() in metabolic_targets:
                node.set_fillcolor('#8EC127')

            elif node.get_name() in sites_fc:
                node.set_fillcolor('#3498db')

            elif node.get_name() in kinases_es:
                node.set_fillcolor('#BB3011')

            node.set_name(uniprot2gene[node.get_name()][0] if node.get_name() in uniprot2gene else node.get_name())
            graph.add_node(node)

        edge = pydot.Edge(source, target)
        graph.add_edge(edge)

graph.write_pdf('/Users/emanuel/Downloads/test.pdf')
print '[INFO] Network PDF saved!'

# ABL1: P00519, PDK1: Q15118, PDHA1: P08559, PDP1: Q9P0J1, ABL2: P42684, LASP1: Q14847