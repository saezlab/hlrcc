import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from hlrcc.utils.volcano import volcano
from pymist.reader.sbml_reader import read_sbml_model
from pymist.utils.map_peptide_sequence import read_fasta, read_uniprot_genename
from pandas import DataFrame, Series, read_csv


# -- Import data-sets
tp = read_csv('%s/data/uok262_proteomics_logfc.txt' % wd, sep='\t')
pp = read_csv('%s/data/uok262_phosphoproteomics_logfc.txt' % wd, sep='\t')


# -- Label uniprot ids
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)

tp['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else '' for i in tp.index]
pp['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else '' for i in pp.index]


# -- Import metabolic model
m_genes = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml').get_genes()


# -- Plot volcanos
genes_highlight = ['VIM', 'PDHA1', 'GAPDH', 'FH', 'ABL1', 'ABL2']

volcano(
    '%s/reports/phosphoproteomics_logfc_volcano.pdf' % wd,
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
    '%s/reports/proteomics_logfc_volcano.pdf' % wd,
    tp,
    'logFC',
    'p.value.log10',
    'adj.P.Val',
    'Human - UOK 262 vs UOK 262 pFH \n proteomics',
    genes_highlight,
    'name'
)
plt.close('all')