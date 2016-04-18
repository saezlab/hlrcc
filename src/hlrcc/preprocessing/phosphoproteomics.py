import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv
from hlrcc import data_dir, wd
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
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'pp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index


# -- Preprocess
# Import and process phospho
info_columns = ['peptide', 'site', 'uniprot']

pp_all = read_csv('%s/data/uok262_phosphoproteomics.txt' % wd, sep='\t').dropna(subset=info_columns)
pp = pp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', pp.shape

# Remove peptides with 1 or less measurements per condition
pp = pp[(pp[ss_ko].count(1) > 1) & (pp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', pp.shape

# Considering proteotypic peptides
pp = pp[[len(i.split('; ')) == 2 for i in pp['uniprot']]]
print '[INFO] Considering proteotypic peptides: ', pp.shape

# Consider single phosphorylated peptides
pp = pp[[len(i.split('+')) == 1 for i in pp['site']]]
print '[INFO] Consider singly phosphorylated peptides: ', pp.shape

# Create p-site IDs
pp['index'] = pp.index
pp['psite'] = ['%s_%s_%s_%d' % (prot.split('; ')[0], re.findall('\(([A-Z][0-9]*)\)', psite)[0], pep, index) for prot, psite, pep, index in pp[['uniprot', 'site', 'peptide', 'index']].values]
pp = pp.drop(info_columns + ['index'], axis=1).set_index('psite')

# Log 2 transform
pp[ss.index] = np.log2(pp[ss.index])

# Scale samples
pp[ss.index] = zscore(pp[ss.index])

# Export p-site level phosphoproteomics
pp.columns = [ss.ix[i, 'condition'] for i in pp.columns]
pp.to_csv('%s/data/uok262_phosphoproteomics_processed.txt' % wd, sep='\t')
print '[INFO] Export p-site level phosphoproteomics'


# -- Plot
# Heatmap
cmap = sns.light_palette((210, 90, 60), input='husl', as_cmap=True)

sns.clustermap(pp.corr(method='pearson'), annot=True, cmap=cmap)
plt.savefig('%s/reports/phosphoproteomics_replicates_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'

# Volcano
pp_fc = read_csv('%s/data/uok262_phosphoproteomics_logfc.txt' % wd, sep='\t')
pp_fc['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else '' for i in pp_fc.index]

genes_highlight = ['VIM', 'PDHA1', 'GAPDH', 'FH', 'ABL1', 'ABL2']

volcano(
    '%s/reports/phosphoproteomics_logfc_volcano.pdf' % wd,
    pp_fc,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'UOK262 phosphoproteomics (KO vs WT)',
    genes_highlight,
    'name'
)
plt.close('all')
