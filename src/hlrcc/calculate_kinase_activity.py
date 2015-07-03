import pickle
import numpy as np
from hlrcc import wd
from pymist.enrichment.gsea import gsea
from pandas import read_csv, Series

# Import phospho data
pp_fc = read_csv(wd + '/data/b1368p100_phospho_human_limma.tsv', sep='\t', index_col=0)['logFC']

# ---- Calculate kinase activity scores
# Import network
ks_network = read_csv('/Users/emanuel/Projects/resources/signalling_models/KS_network.tab', sep='\t', header=0, dtype=str).dropna()
ks_network['RID'] = [x[0] for x in ks_network['SID']]
ks_network['site'] = ks_network['res'] + ks_network['pos']
print '[INFO] network: ', ks_network.shape

# Filter predicted interactions
ks_network = ks_network.loc[[not x.startswith('n') for x in ks_network['SID']]]
print '[INFO] network, remove predicted interactions: ', ks_network.shape

# Remove self phosphorylations
ks_network = ks_network[[r['S.AC'] != r['K.AC'] for i, r in ks_network.iterrows()]]
print '[INFO] network, remove self phosphorylations: ', ks_network.shape

# Calculate kinase target sites
kinases = set(ks_network.loc[ks_network['K.AC'] != '', 'K.AC'])
kinases_targets = {k: set(map('_'.join, ks_network.loc[ks_network['K.AC'] == k, ['S.AC', 'site']].values)) for k in kinases}

with open('%s/files/kinases_targets.pickle' % wd, 'wb') as handle:
    pickle.dump(kinases_targets, handle)

kinases_targets = {k: kinases_targets[k].intersection(pp_fc.index) for k in kinases_targets}
kinases_targets = {k: v for k, v in kinases_targets.iteritems() if len(v) > 0}
print '[INFO] Kinases targets: ', len(kinases_targets)

# Calculate kinase activity score
kinases_es = {k: gsea(pp_fc.to_dict(), kinases_targets[k], 10000) for k in kinases_targets}
kinases_es = {k: -np.log10(kinases_es[k][1]) if kinases_es[k][0] < 0 else np.log10(kinases_es[k][1]) for k in kinases_es}
Series(kinases_es).to_csv(wd + '/files/kinase_activity.tab', sep='\t')
print '[INFO] Kinase enrichment: ', len(kinases_es)