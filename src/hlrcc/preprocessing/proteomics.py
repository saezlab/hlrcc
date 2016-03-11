import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pandas import read_csv
from pandas.stats.misc import zscore
from hlrcc.utils.volcano import volcano
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.reader.sbml_reader import read_sbml_model


# -- Uniprot names
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)


# -- Import metabolic model
m_genes = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml').get_genes()


# -- Import samplesheet
ss = read_csv('%s/data/proteomics_samplesheet.txt' % wd, sep='\t', index_col=0)
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'tp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index


# -- Process
# Import and process phospho
info_columns = ['peptide', 'uniprot']

tp_all = read_csv('%s/data/uok262_proteomics.txt' % wd, sep='\t').dropna(subset=info_columns)
tp = tp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', tp.shape

# Remove peptides with 1 or less measurements per condition
tp = tp[np.bitwise_and(tp[ss_ko].count(1) > 1, tp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', tp.shape

# Considering proteotypic peptides
tp = tp[[(len(i.split('; ')) == 2) and i.split('; ')[0] != '' for i in tp['uniprot']]]
print '[INFO] Considering proteotypic peptides: ', tp.shape

# Create protein IDs
tp['uniprot'] = [i.split(';')[0] for i in tp['uniprot']]
tp = tp.groupby('uniprot', sort=False).median()

# Log 2 transform
tp[ss.index] = np.log2(tp[ss.index])

# Scale samples
tp[ss.index] = zscore(tp[ss.index])

# Export protein level proteomics
tp.columns = [ss.ix[i, 'condition'] for i in tp.columns]
tp.to_csv('%s/data/uok262_proteomics_processed.txt' % wd, sep='\t')
print '[INFO] Export protein level proteomics'


# -- Plot
# Heatmap
cmap = sns.light_palette((210, 90, 60), input='husl', as_cmap=True)

sns.clustermap(tp.corr(method='pearson'), annot=True, cmap=cmap)
plt.savefig('%s/reports/proteomics_replicates_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'

# Volcano
tp_fc = read_csv('%s/data/uok262_proteomics_logfc.txt' % wd, sep='\t')
tp_fc['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else '' for i in tp_fc.index]

genes_highlight = ['VIM', 'PDHA1', 'GAPDH', 'FH', 'ABL1', 'ABL2']

volcano(
    '%s/reports/proteomics_logfc_volcano.pdf' % wd,
    tp_fc,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'UOK262 proteomics (KO vs WT)',
    genes_highlight,
    'name'
)
plt.close('all')
