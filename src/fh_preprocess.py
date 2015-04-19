__author__ = 'emanuel'

import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
from pandas import read_csv, DataFrame
from utils.map_peptide_sequence import read_fasta, match_sequence, read_uniprot_accname, read_uniprot_genename
from sklearn.decomposition import PCA

# Configure vars
data_dir = '/Users/emanuel/Projects/data/fh_cells/'
organism = ['human', 'mouse'][0]

os = {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[organism]
acc2uniprot = read_uniprot_accname(os=os)
uniprot_fasta = read_fasta(os=os)
uniprot2genename = read_uniprot_genename(os=os)

# Define files paths
pp_file = {
    'human': 'human_phosphoproteomics/b1368p100_phosho_human.tab',
    'mouse': 'mouse_phosphoproteomics/b1368p100_phosho_mouse.tab'
}[organism]

tp_file = {
    'human': 'human_proteomics/b1368p100_protein_human.tab',
    'mouse': 'mouse_proteomics/b1368p100_protein_mouse.tab'
}[organism]

pp_file_processed = {
    'human': 'human_phosphoproteomics/b1368p100_phosho_human_processed.tab',
    'mouse': 'mouse_phosphoproteomics/b1368p100_phosho_mouse_processed.tab'
}[organism]

tp_file_processed = {
    'human': 'human_proteomics/b1368p100_protein_human_processed.tab',
    'mouse': 'mouse_proteomics/b1368p100_protein_mouse_processed.tab'
}[organism]

# Import samplesheet
ss = read_csv(data_dir + 'fh_samplesheet.tab', sep='\t', index_col=0)
ss = ss.loc[ss['organism'] == organism]

ss_pp = list(ss[ss['type'] == 'pp'].index)
ss_pp_ko = list(ss[(ss['type'] == 'pp') & (ss['condition'] == 'fh_ko')].index)
ss_pp_wt = list(ss[(ss['type'] == 'pp') & (ss['condition'] == 'fh_wt')].index)

ss_tp = list(ss[ss['type'] == 'tp'].index)
ss_tp_ko = list(ss[(ss['type'] == 'tp') & (ss['condition'] == 'fh_ko')].index)
ss_tp_wt = list(ss[(ss['type'] == 'tp') & (ss['condition'] == 'fh_wt')].index)

#### Import and preprocess phospho and proteomics data
pp_all = read_csv(data_dir + pp_file, sep='\t')
tp_all = read_csv(data_dir + tp_file, sep='\t')

# Drop NaN on phospho sites and peptides
pp_all = pp_all.dropna(subset=['peptide', 'site'])
tp_all = tp_all.dropna(subset=['peptide'])

# Assemble reduced data-set
pp_columns = ['peptide', 'site', 'uniprot']
pp_columns.extend(ss_pp)
pp = pp_all[pp_columns]

tp_columns = ['peptide', 'acc_no', 'uniprot']
tp_columns.extend(ss_tp)
tp = tp_all[tp_columns]

# If mouse data rematch peptide with mouse uniprot Ids
if organism == 'mouse':
    pp['uniprot'] = [match_sequence(uniprot_fasta, x) for x in pp['peptide']]
    pp['uniprot'] = ['; '.join(x) for x in pp['uniprot']]

# Replace zeros with NaN
pp = pp.replace(0.0, np.NaN)
tp = tp.replace(0.0, np.NaN)

# Remove peptides with 1 or less measurements per condition
pp = pp.loc[(np.isnan(pp.ix[:, ss_pp_ko]).sum(axis=1) < 5) & (np.isnan(pp.ix[:, ss_pp_wt]).sum(axis=1) < 5), ]
tp = tp.loc[(np.isnan(tp.ix[:, ss_tp_ko]).sum(axis=1) < 2) & (np.isnan(tp.ix[:, ss_tp_wt]).sum(axis=1) < 2), ]

#### Log 2 transform
tp[ss_tp] = np.log2(tp[ss_tp])
pp[ss_pp] = np.log2(pp[ss_pp])

#### Scale replicates
tp[ss_tp] = (tp[ss_tp] - tp[ss_tp].mean()) / tp[ss_tp].std()
pp[ss_pp] = (pp[ss_pp] - pp[ss_pp].mean()) / pp[ss_pp].std()

#### Export
pp.to_csv(data_dir + pp_file_processed, sep='\t', index=False)
tp.to_csv(data_dir + tp_file_processed, sep='\t', index=False)

#### Verbose
print '[INFO] Preprocess done!'