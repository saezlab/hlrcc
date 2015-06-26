import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from hlrcc.utils.volcano import volcano
from pymist.reader.sbml_reader import read_sbml_model
from bioservices import UniProt
from pandas import DataFrame, Series, read_csv


def get_fasta(i):
    try:
        return re.findall('.* GN=(.*?) ', uniprot_service.get_fasta(i))

    except TypeError:
        return []

# ---- Import data-sets
tp = read_csv('%s/data/b1368p100_protein_human_limma.tsv' % wd, sep='\t')
pp = read_csv('%s/data/b1368p100_phospho_human_limma.tsv' % wd, sep='\t')

# ---- Import uniprot information
uniprot_service = UniProt(cache=True)

uniprot_map = set(np.concatenate((tp.index, [i.split('_')[0] for i in pp.index])))
uniprot_map = {i: get_fasta(i) for i in uniprot_map}
uniprot_map = {i: uniprot_map[i][0] if len(uniprot_map[i]) == 1 else '' for i in uniprot_map}
print '[INFO] UniProt ID and Gene Name map loaded'

tp['name'] = [uniprot_map[i] for i in tp.index]

pp['name'] = [i.split('_')[0] for i in pp.index]
pp['name'] = [uniprot_map[i] for i in pp['name']]

# ---- Import metabolic model
model = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml')
m_genes = model.get_genes()

# ---- Plot volcanos
genes_highlight = ['VIM', 'PDHA1', 'GAPDH', 'FH', 'ABL1']

volcano(
    '%s/reports/phosphoproteomics_volcano.pdf' % wd,
    pp,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - UOK 262 vs UOK 262 pFH \n phosphoproteomics',
    genes_highlight,
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
    genes_highlight,
    'name'
)
plt.close('all')