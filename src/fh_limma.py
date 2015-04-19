__author__ = 'emanuel'

import sys
sys.path.extend(['/Users/emanuel/Projects/projects/pymist', '/Users/emanuel/Projects/projects/pymist/pymist'])

import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv
from utils.map_peptide_sequence import read_uniprot_genename, read_fasta, match_sequence


def cohensd(c0, c1):
    return (np.mean(c0) - np.mean(c1)) / (np.sqrt((np.std(c0) ** 2 + np.std(c1) ** 2) / 2))


def signif(value):
    if value < 0.01:
        return '*'
    elif value < 0.05:
        return '**'
    else:
        return '-'


def plot_volcanos_results(f, d, title='', genes_to_highlight=None, w=6, h=9):
    sns.set(style='ticks')
    sns.set_context('paper')
    sns.despine(offset=10)

    # Define pallete
    colour_pallete = [sns.color_palette('OrRd')[4], sns.color_palette('OrRd')[3], sns.color_palette('OrRd')[1]]
    sns.lmplot(x='logFC', y='p.value.log10', data=d, hue='signif', fit_reg=False, palette=colour_pallete, legend=False)

    # Add FDR threshold lines
    plt.text(plt.xlim()[0]*.98, -np.log10(0.01), 'FDR 1%', ha='left', color=colour_pallete[0], alpha=0.65, fontsize=7)
    plt.axhline(-np.log10(0.01), c=colour_pallete[0], ls='--', lw=.5, alpha=.7)

    plt.text(plt.xlim()[0]*.98, -np.log10(0.05), 'FDR 5%', ha='left', color=colour_pallete[1], alpha=0.65, fontsize=6)
    plt.axhline(-np.log10(0.05), c=colour_pallete[1], ls='--', lw=.5, alpha=.7)

    # Add axis lines
    plt.axvline(0, c='#95a5a6', lw=.3, alpha=.15)
    [plt.axhline(x, c='#95a5a6', lw=.3, alpha=.15) for x in range(0, int(plt.ylim()[1]), 2)]

    # Add axis labels and title
    plt.title(title, fontsize=10, fontname='Arial')
    plt.xlabel('fold change (log)', fontsize=8, fontname='Arial')
    plt.ylabel('p-value (-log)', fontsize=8, fontname='Arial')

    # Add text to highlighted genes
    if genes_to_highlight is not None:
        for i, r in d.iterrows():
            gene = r['gene_symbol']
            if len(set(genes_to_highlight).intersection(set(gene))) >= 1:
                plt.text(r['logFC'] * 1.01, r['p.value.log10'] * 1.01, '; '.join(gene), ha='left', alpha=0.75, fontsize=8)

    # Adjust axis lines thickness
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.xaxis.set_tick_params(width=0.5)
    ax.yaxis.set_tick_params(width=0.5)

    # Save plot
    fig = plt.gcf()
    fig.set_size_inches(w, h)
    fig.savefig(f)

    # Clean plot
    plt.close()

    print '[INFO] Volcano generated: ' + f

# Configure vars
results_dir = '/Users/emanuel/Projects/projects/pymist/reports/fh_cells/'
uniprot_dir = '/Users/emanuel/Projects/resources/uniprot/'

organism = ['human', 'mouse'][0]
os = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[organism]
uniprot2gene = read_uniprot_genename(os=os)
uniprot_seq = read_fasta(os=os)

# Define files paths
pp_results_file = results_dir + {
    'human': 'b1368p100_phospho_human_limma.tsv',
    'mouse': 'b1368p100_phospho_mouse_limma.tsv'
}[organism]

tp_results_file = results_dir + {
    'human': 'b1368p100_protein_human_limma.tsv',
    'mouse': 'b1368p100_protein_mouse_limma.tsv'
}[organism]

pp_volcano_file = results_dir + {
    'human': time.strftime('%X') + '_' + 'human_phospho_volcano',
    'mouse': time.strftime('%X') + '_' + 'mouse_phospho_volcano',
}[organism]

tp_volcano_file = results_dir + {
    'human': time.strftime('%X') + '_' + 'human_protein_volcano',
    'mouse': time.strftime('%X') + '_' + 'mouse_protein_volcano',
}[organism]

# Import limma differential analysis
pp_results = read_csv(pp_results_file, sep='\t')
tp_results = read_csv(tp_results_file, sep='\t')

# Descritise significance
pp_results['signif'] = [signif(v) for v in pp_results['adj.P.Val']]
tp_results['signif'] = [signif(v) for v in tp_results['adj.P.Val']]

# Convert uniprot to gene symbol
# pp_results['gene_symbol'] = [uniprot2gene[x] if x in uniprot2gene else '' for x in pp_results['uniprot']]
# tp_results['gene_symbol'] = [uniprot2gene[x] if x in uniprot2gene else '' for x in tp_results['uniprot']]

# Map peptide to uniprot
pp_results['uniprots'] = [match_sequence(uniprot_seq, pep) for pep in pp_results['peptide']]
tp_results['uniprots'] = [match_sequence(uniprot_seq, pep) for pep in tp_results['peptide']]

# Remove very ambiguos peptides
pp_results = pp_results[[len(x) < 4 for x in pp_results['uniprots']]]
tp_results = tp_results[[len(x) < 4 for x in tp_results['uniprots']]]

# Map uniprot to gene symbol
pp_results['gene_symbol'] = [[uniprot2gene[x] for x in u if x in uniprot2gene] for u in pp_results['uniprots']]
tp_results['gene_symbol'] = [[uniprot2gene[x] for x in u if x in uniprot2gene] for u in tp_results['uniprots']]

# Plot volcanos
genes_to_highlight = {
    'human': {'VIM', 'FH', 'PDHA1', 'GAPDH', 'LASP1'},
    'mouse': {'Vim', 'Fh', 'Rai14', 'Ckb', 'Pdha1', 'Gapdh', 'Des', 'Gsta4', 'Idh1', 'Ablim1', 'Eif4b'},
}[organism]

plot_volcanos_results(pp_volcano_file + '.pdf', pp_results, 'Phosphoproteomics - KO vs WT (' + organism + ')', genes_to_highlight)
plot_volcanos_results(tp_volcano_file + '.pdf', tp_results, 'Proteomics - KO vs WT (' + organism + ')', genes_to_highlight)

# # Protein behaviour
# pp_uniprot = [uniprot for uniprots in pp_results['uniprots'] for uniprot in uniprots if len(uniprots) < 2]
# tp_uniprot = [uniprot for uniprots in tp_results['uniprots'] for uniprot in uniprots]
#
# overlap_uniprot = list(set(pp_uniprot).intersection(set(tp_uniprot)))
#
# stats, fc_thres = {}, 0.6
# for uniprot in overlap_uniprot:
#     pp_values = pp_results[[uniprot in values for values in pp_results['uniprots']]]
#     tp_values = tp_results[[uniprot in values for values in tp_results['uniprots']]]
#
#     pp_values_median = np.median(pp_values['logFC'])
#     tp_values_median = np.median(tp_values['logFC'])
#
#     if abs(tp_values_median) < fc_thres and abs(pp_values_median) < fc_thres:
#         stats[uniprot] = 1
#
#     elif (tp_values_median > fc_thres and pp_values_median > fc_thres) or (tp_values_median < -fc_thres and pp_values_median < -fc_thres):
#         stats[uniprot] = 2
#
#     elif abs(tp_values_median) > fc_thres and abs(pp_values_median) < fc_thres:
#         stats[uniprot] = 3
#
#     elif abs(tp_values_median) < fc_thres and abs(pp_values_median) > fc_thres:
#         stats[uniprot] = 4
#
#     elif (tp_values_median > fc_thres and pp_values_median < -fc_thres) or (tp_values_median < -fc_thres and pp_values_median > fc_thres):
#         stats[uniprot] = 5
#
# behaviour_type = Series(stats.values()).value_counts()
# behaviour_type.plot(kind='pie'); plt.show();