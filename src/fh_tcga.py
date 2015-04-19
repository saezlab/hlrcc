__author__ = 'emanuel'

import scipy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv


types = ['KICH', 'KIRC', 'KIRP']

wd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/'

genes_to_highlight = {'FH', 'VIM', 'CDH1'}

for tumour in types:
    tumour_df = read_csv(wd + 'files/tcga_' + tumour + '_limma_rsem_raw_counts.tab', sep='\t')
    tumour_df['gene'] = tumour_df.index

    fig = plt.figure()

    sns.set_style('white')
    sns.barplot('gene', 'logFC', data=tumour_df.ix[genes_to_highlight])
    sns.despine()

    fig.set_size_inches(5, 7)

    plt.axhline(0, c='#95a5a6', lw=.3, alpha=.6)
    plt.title('Tumour vs Normal - ' + tumour)

    plt.savefig(wd + 'plots/tcga_' + tumour + '_barplot.pdf')
    plt.close('all')

    tumour_df = read_csv('/Users/emanuel/Projects/projects/tcga_rna_seq/data_preprocessed/' + tumour + '_rsem_tumour.tab', sep='\t', index_col=0)

    sns.jointplot(tumour_df.ix['FH'], tumour_df.ix['CDH1'], kind='reg', size=5, ratio=6)
    plt.savefig(wd + 'plots/tcga_' + tumour + '_FH_CDH1_cor.pdf', bbox_inches='tight')
    plt.close('all')

    sns.jointplot(tumour_df.ix['FH'], tumour_df.ix['VIM'], kind='reg', size=5, ratio=6)
    plt.savefig(wd + 'plots/tcga_' + tumour + '_FH_VIM_cor.pdf', bbox_inches='tight')
    plt.close('all')