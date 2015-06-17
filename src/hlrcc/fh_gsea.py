__author__ = 'emanuel'

import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename
from pymist.enrichment.gsea import gsea, plot_gsea

organism = 1  # 0: human, 1: mouse

wd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/'

# Import data-sets
trans = read_csv([
    wd + 'files/human_transcriptomics.tab',
    wd + 'files/mouse_transcriptomics.tab'
][organism], sep='\t')

# Import signatures
emt = read_csv([
  wd + 'files/human_byers_emt.tab',
  wd + 'files/mouse_byers_emt.tab'
][organism], sep='\t')

for signature in set(emt['signature']):
    escore, pvalue, info = gsea(trans['logFC'], set(emt.loc[emt['signature'] == signature, 'gene']), True, 1000)
    plot_gsea(escore, pvalue, info, wd + 'plots/', 'emt_' + ['human', 'mouse'][organism] + '.pdf', signature)