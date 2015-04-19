__author__ = 'emanuel'

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame
from pymist.enrichment.gsea import gsea
from pymist.utils.map_peptide_sequence import read_uniprot_accname

# Import lists
human_gl = read_csv('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/' + 'etc_proteins_human.tab', sep='\t')

# Import data-sets
human_tp = read_csv('/Users/emanuel/Projects/data/fh_cells/petros_tmt/' + 'human_proteomics.tab', sep='\t', index_col='Accession')

title = 'human_proteomics'

for signature in set(human_gl['signature']):
    signature_genes = set(human_gl.loc[human_gl['signature'] == signature, 'accession'])

    escore, pvalue, info = gsea(human_tp.mean(1).to_dict(), signature_genes, True, 10000)

    x, y = np.array(range(len(info[0]))), np.array(info[0])

    sns.set_style('white')
    plt.plot(x, y, '-')
    sns.despine()
    [plt.axvline(hit, c='#f7941e', lw=.3, alpha=.45) for hit in x[info[1]]]
    plt.axhline(0.0, c='#95a5a6', lw=.3, alpha=.85)
    plt.ylabel('Enrichment score')
    plt.xlim([0, len(info[0])])
    plt.annotate('p-value: ' + str(pvalue), xy=(0.99, 0.95), xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')
    plt.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/petros_tmt/' + title + '_' + signature + '_gsea.pdf', bbox_inches='tight')
    plt.title(title)
    plt.close('all')
    print '[INFO] ' + signature


# Mouse
uniprot_ac = dict((v, k) for k, v in read_uniprot_accname(os='Mus musculus').iteritems())

mouse_tp = read_csv('/Users/emanuel/Projects/data/fh_cells/petros_tmt/' + 'mouse_proteomics.tab', sep='\t', index_col='acc_name').dropna()
mouse_gl = read_csv('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/' + 'etc_proteins_mouse.tab', sep='\t')

mouse_tp['uniprot'] = [uniprot_ac[i] if i in uniprot_ac else np.NaN for i in mouse_tp.index]
mouse_tp = mouse_tp.groupby('uniprot').median()
mouse_tp.to_csv('/Users/emanuel/Projects/data/fh_cells/petros_tmt/' + 'mouse_proteomics_v2.tab', sep='\t')

mouse_gl['uniprot'] = [uniprot_ac[i] for i in mouse_gl['acc_name']]
mouse_gl.to_csv('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/' + 'etc_proteins_mouse_v2.tab', sep='\t')

title = 'mouse_proteomics'

for signature in set(mouse_gl['signatue']):
    signature_genes = set(mouse_gl.loc[mouse_gl['signatue'] == signature, 'uniprot'])

    escore, pvalue, info = gsea(mouse_tp.mean(1).to_dict(), signature_genes, True, 10000)

    x, y = np.array(range(len(info[0]))), np.array(info[0])

    sns.set_style('white')
    plt.plot(x, y, '-')
    sns.despine()
    [plt.axvline(hit, c='#f7941e', lw=.3, alpha=.45) for hit in x[info[1]]]
    plt.axhline(0.0, c='#95a5a6', lw=.3, alpha=.85)
    plt.ylabel('Enrichment score')
    plt.xlim([0, len(info[0])])
    plt.annotate('p-value: ' + str(pvalue), xy=(0.99, 0.95), xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')
    plt.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/petros_tmt/' + title + '_' + signature + '_gsea.pdf', bbox_inches='tight')
    plt.title(title)
    plt.close('all')
    print '[INFO] ' + signature

