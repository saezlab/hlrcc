import re
import numpy as np
from pandas.stats.misc import zscore
from hlrcc import data_dir, wd
from pandas import read_csv

# ---- Import samplesheet
ss = read_csv(data_dir + '/fh_samplesheet.tab', sep='\t', index_col=0)
ss = ss.loc[np.bitwise_and(ss['organism'] == 'human', ss['type'] == 'pp')]

ss_ko = ss[ss['condition'] == 'fh_ko'].index
ss_wt = ss[ss['condition'] == 'fh_wt'].index

# ---- Process
# Import and process phospho
info_columns = ['peptide', 'site', 'uniprot']

pp_all = read_csv(wd + '/data/b1368p100_phosho_human.tab', sep='\t').dropna(subset=info_columns)
pp = pp_all[np.concatenate((info_columns, ss.index))].replace(0.0, np.NaN)
print '[INFO] phospho: ', pp.shape

# Remove peptides with 1 or less measurements per condition
pp = pp[np.bitwise_and(pp[ss_ko].count(1) > 1, pp[ss_wt].count(1) > 1)]
print '[INFO] Remove peptides with 1 or less measurements per condition: ', pp.shape

# Log 2 transform
pp[ss.index] = np.log2(pp[ss.index])

# Scale samples
pp[ss.index] = zscore(pp[ss.index])

# Considering proteotypic peptides
pp = pp[[len(i.split('; ')) == 2 for i in pp['uniprot']]]
pp['uniprot'] = [i.split('; ')[0] for i in pp['uniprot']]
print '[INFO] Considering proteotypic peptides: ', pp.shape

# Consider singly phosphorylated peptides
pp = pp[[len(i.split('+')) == 1 for i in pp['site']]]
pp['site'] = [re.findall('\(([A-Z][0-9]*)\)', i)[0] for i in pp['site']]
print '[INFO] Consider singly phosphorylated peptides: ', pp.shape

# Average p-site phosphorylation
pp['site'] = pp['uniprot'] + '_' + pp['site']
pp = pp.groupby('site').median()
print '[INFO] Average p-site phosphorylation: ', pp.shape

# Export p-site level phosphoproteomics
pp.columns = [ss.ix[i, 'condition'] for i in pp.columns]
pp.to_csv(wd + '/data/b1368p100_phospho_human_processed.tsv', sep='\t')
print '[INFO] Export p-site level phosphoproteomics'