import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from hlrcc.volcano import volcano
from bioservices import UniProt
from pandas import DataFrame, Series, read_csv

# ---- Import data-sets
tp = read_csv('%s/files/b1368p100_protein_human_limma.tsv' % wd, sep='\t')
pp = read_csv('%s/files/b1368p100_phospho_human_limma.tsv' % wd, sep='\t')

# ---- Import uniprot information
uniprot_service = UniProt(cache=True)
uniprot_info = {i: uniprot_service.get_df(i) for i in set(np.concatenate((tp.index, [i.split('_')[0] for i in pp.index])))}
print '[INFO] Bioservices loaded unirpto info: ', len(uniprot_info)

# ---- Plot volcanos
volcano(
    '%s/reports/phosphoproteomics_volcano.pdf' % wd,
    pp,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - UOK 262 vs UOK 262 pFH \n phosphoproteomics'
)