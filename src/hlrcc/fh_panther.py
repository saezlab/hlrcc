__author__ = 'emanuel'

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv, DataFrame
from utils.map_peptide_sequence import read_uniprot_accname, read_uniprot_genename
from utils.plot_utils import save_plot

organism = ['human', 'mouse'][0]

mito_complexes = [
    'Mitochondrial_complex_I.txt',
    'Mitochondrial_complex_II.txt',
    'Mitochondrial_complex_III.txt',
    'Mitochondrial_complex_IV.txt',
    'Mitochondrial_complex_V.txt'
]

# Configure vars
results_dir = '/Users/emanuel/Projects/projects/pymist/reports/fh_cells/'
uniprot_dir = '/Users/emanuel/Projects/resources/uniprot/'

os = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[organism]
acc2uniprot = read_uniprot_accname(os=os)
acc2uniprot = dict(zip(acc2uniprot.values(), acc2uniprot.keys()))

uniprot2gene = read_uniprot_genename(os=os)

# Define files paths
pp_results_file = results_dir + {
    'human': 'b1368p100_phospho_human_limma.tsv',
    'mouse': 'b1368p100_phospho_mouse_limma.tsv'
}[organism]

tp_results_file = results_dir + {
    'human': 'b1368p100_protein_human_limma.tsv',
    'mouse': 'b1368p100_protein_mouse_limma.tsv'
}[organism]

# Import limma differential analysis
pp_results = read_csv(pp_results_file, sep='\t')
tp_results = read_csv(tp_results_file, sep='\t')

for pathway_file in mito_complexes:

    pathway = read_csv('/Users/emanuel/Projects/resources/panther/' + pathway_file, sep='\t')
    pathway_genes = [g.upper() for g in pathway['Approved Symbol']]

    tp_pathway = [(uniprot2gene[acc2uniprot[v['acc_no']]], v['logFC']) for gene in pathway_genes for k, v in tp_results.iterrows() if v['acc_no'] in acc2uniprot and acc2uniprot[v['acc_no']] in uniprot2gene and uniprot2gene[acc2uniprot[v['acc_no']]].upper() == gene.upper()]
    tp_pathway = DataFrame(tp_pathway, columns=['uniprot', 'logFC'])

    # # Plot pathway measurements
    # panther_pathway_uniprot = {v[4].split('UniProtKB=')[1]: v[3] for k, v in panther_pathway.iterrows()}
    # tp_pathway = [(uniprot2gene[uniprot], v['logFC']) for uniprot in panther_pathway_uniprot for k, v in tp_results.iterrows() if v['acc_no'] in acc2uniprot and acc2uniprot[v['acc_no']] == uniprot]
    # tp_pathway = DataFrame(tp_pathway, columns=['uniprot', 'logFC'])

    sns.set(style='white')
    g = sns.factorplot('uniprot', 'logFC', data=tp_pathway, hue='uniprot', kind='point')
    plt.axhline(0, c='#95a5a6', lw=.3, alpha=.85)
    plt.title(organism.upper() + ': KO vs WT')
    g.set_xticklabels(rotation=30)
    plt.ylabel('Peptide intensities logFC')
    save_plot('/Users/emanuel/Downloads/' + organism + '_' + 'total_protein.pdf', 10, 5)