import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pandas import read_csv
from pandas.stats.misc import zscore


# -- Import samplesheet
ss = read_csv('%s/data/fh_samplesheet.tab' % wd, sep='\t', index_col=0)
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'tp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index


# -- Process
# Import and process phospho
info_columns = ['peptide', 'uniprot']

tp_all = read_csv('%s/data/uok262_proteomics.tab' % wd, sep='\t').dropna(subset=info_columns)
tp = tp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', tp.shape

# Remove peptides with 1 or less measurements per condition
tp = tp[np.bitwise_and(tp[ss_ko].count(1) > 1, tp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', tp.shape

# Considering proteotypic peptides
tp = tp[[(len(i.split('; ')) == 2) and i.split('; ')[0] != '' for i in tp['uniprot']]]
print '[INFO] Considering proteotypic peptides: ', tp.shape

# Create protein IDs
tp['index'] = tp.index
tp['protein'] = ['%s_%s_%d' % (uniprot.split(';')[0], pep, index) for pep, uniprot, index in tp[['peptide', 'uniprot', 'index']].values]
tp = tp.drop(info_columns + ['index'], axis=1).set_index('protein')

# Log 2 transform
tp[ss.index] = np.log2(tp[ss.index])

# Scale samples
tp[ss.index] = zscore(tp[ss.index])

# Export protein level proteomics
tp.columns = [ss.ix[i, 'condition'] for i in tp.columns]
tp.to_csv('%s/data/uok262_proteomics_processed.tab' % wd, sep='\t')
print '[INFO] Export protein level proteomics'


# -- Plot heatmap
cmap = sns.light_palette((210, 90, 60), input='husl', as_cmap=True)

sns.clustermap(tp.corr(method='pearson'), annot=True, cmap=cmap)
plt.savefig('%s/reports/proteomics_replicates_clustermap.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Heatmap plotted!'
