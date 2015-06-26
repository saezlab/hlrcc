import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from hlrcc.utils.volcano import volcano
from bioservices import UniProt
from pandas import DataFrame, Series, read_csv

# ---- Import data-sets
tp = read_csv('%s/files/b1368p100_protein_human_limma.tsv' % wd, sep='\t')
pp = read_csv('%s/files/b1368p100_phospho_human_limma.tsv' % wd, sep='\t')

# ---- Import uniprot information
uniprot_service = UniProt(cache=True)
uniprot_info = {}
for i in set(np.concatenate((tp.index, [i.split('_')[0] for i in pp.index]))):
    try:
        res = uniprot_service.get_df(i, organism='human')
        uniprot_info[i] = res

    except AttributeError:
        uniprot_info[i] = ''
        print '[Warning] Uniprot ID %s not found' % i

print '[INFO] Bioservices loaded unirpto info: ', len(uniprot_info)

tp['name'] = [uniprot_info[i]['Gene names'][0][0].split(' ')[0] if i != '' and len(uniprot_info[i]['Gene names'][0]) > 0 else '' for i in tp.index]

pp['name'] = [i.split('_')[0] for i in pp.index]
pp['name'] = [uniprot_info[i]['Gene names'][0][0].split(' ')[0] if i != '' and len(uniprot_info[i]['Gene names'][0]) > 0 else '' for i in pp['name']]

# ---- Plot volcanos
volcano(
    '%s/reports/phosphoproteomics_volcano.pdf' % wd,
    pp,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - UOK 262 vs UOK 262 pFH \n phosphoproteomics',
    ['VIM', 'PDHA1', 'GAPDH', 'FH', 'LASP1'],
    'name'
)
plt.close('all')

volcano(
    '%s/reports/proteomics_volcano.pdf' % wd,
    tp,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - UOK 262 vs UOK 262 pFH \n proteomics',
    ['VIM', 'PDHA1', 'GAPDH', 'FH', 'LASP1'],
    'name'
)
plt.close('all')