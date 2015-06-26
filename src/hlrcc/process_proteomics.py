import numpy as np
from pandas.stats.misc import zscore
from hlrcc import wd
from pandas import read_csv

# ---- Import samplesheet
ss = read_csv(wd + '/data/fh_samplesheet.tab', sep='\t', index_col=0)
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'tp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index

# ---- Process
# Import and process phospho
info_columns = ['peptide', 'uniprot']

tp_all = read_csv(wd + '/data/b1368p100_protein_human.tab', sep='\t').dropna(subset=info_columns)
tp = tp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', tp.shape

# Remove peptides with 1 or less measurements per condition
tp = tp[np.bitwise_and(tp[ss_ko].count(1) > 1, tp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', tp.shape

# Log 2 transform
tp[ss.index] = np.log2(tp[ss.index])

# Scale samples
tp[ss.index] = zscore(tp[ss.index])

# Considering proteotypic peptides
tp = tp[[len(i.split('; ')) == 2 for i in tp['uniprot']]]
tp['uniprot'] = [i.split('; ')[0] for i in tp['uniprot']]
tp = tp[[i != '' for i in tp['uniprot']]]
print '[INFO] Considering proteotypic peptides: ', tp.shape

# Average intensities to protein level
tp = tp.groupby('uniprot').median()
print '[INFO] Average intensities to protein level: ', tp.shape

# Export p-site level phosphoproteomics
tp.columns = [ss.ix[i, 'condition'] for i in tp.columns]
tp.to_csv(wd + '/data/b1368p100_protein_human_processed.tab', sep='\t')
print '[INFO] Export protein level proteomics'