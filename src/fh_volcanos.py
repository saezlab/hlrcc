__author__ = 'emanuel'

import sys
sys.path.extend(['/Users/emanuel/Projects/projects/pymist', '/Users/emanuel/Projects/projects/pymist/pymist'])

import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


def signif(value):
    if value < 0.01:
        return '*'
    elif value < 0.05:
        return '**'
    else:
        return '-'


def plot_volcanos_results(f, d, genes_to_highlight, gene_column, title='', w=6, h=9):
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
    for i, r in d.iterrows():
        gene = set([r[gene_column]])
        if len(genes_to_highlight.intersection(gene)) > 0 and r['signif'] != '-':
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

organism = 1  # 0: human, 1: mouse

wd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/'

acc2gene = {acc: gene for uniprot, (gene, acc) in read_uniprot_genename(os=['Homo sapiens', 'Mus musculus'][organism]).items()}

# Import data-sets
trans = read_csv([
    wd + 'files/human_transcriptomics.tab',
    wd + 'files/mouse_transcriptomics.tab'
][organism], sep='\t')
trans['p.value.log10'] = -np.log10(trans['adj.P.Val'])
trans['gene'] = trans.index

prote = read_csv([
    wd + 'files/b1368p100_protein_human_limma.tsv',
    wd + 'files/b1368p100_protein_mouse_limma.tsv'
][organism], sep='\t')
prote['gene'] = [acc2gene[i] if i in acc2gene else '' for i in prote['acc_no']]

volcano_file = [
    'human_volcano.pdf',
    'mouse_volcano.pdf'
][organism]

# Descritise significance
trans['signif'] = [signif(v) for v in trans['adj.P.Val']]
prote['signif'] = [signif(v) for v in prote['adj.P.Val']]

genes_to_highlight = {'Fh1', 'Vim', 'Cdh1'} if organism == 1 else {'FH', 'VIM', 'CDH1'}

plot_volcanos_results(wd + 'transcriptomics_' + volcano_file, trans, genes_to_highlight, 'gene', 'KO vs WT')
plot_volcanos_results(wd + 'proteomics_' + volcano_file, prote, genes_to_highlight, 'gene', 'KO vs WT')

for type in ['KIRC']:
    genes_to_highlight = {'FH', 'VIM', 'CDH1'}

    tcga_rnaseq = read_csv(wd + 'files/tcga_' + type + '_limma_rsem_raw_counts.tab', sep='\t')

    tcga_rnaseq['signif'] = [signif(v) for v in tcga_rnaseq['adj.P.Val']]
    tcga_rnaseq['gene'] = tcga_rnaseq.index
    tcga_rnaseq['p.value.log10'] = -np.log10(tcga_rnaseq['adj.P.Val'])

    plot_volcanos_results(wd + 'plots/tcga_' + type + '_volcano.pdf', tcga_rnaseq, genes_to_highlight, 'gene', 'KO vs WT')