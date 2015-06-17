__author__ = 'emanuel'

import scipy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv


wd, wd2 = '/Users/emanuel/Projects/data/fh_cells/', '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/'

organism = 0  # 0: human, 1: mouse

array_details = [
    wd + 'human_transcriptomics/HuEx-1_0-st-v2.na33.1.hg19.transcript.csv/HuEx-1_0-st-v2.na33.1.hg19.transcript.csv',
    wd + 'mouse_transcriptomics/MoEx-1_0-st-v1.na33.1.mm9.transcript.csv/MoEx-1_0-st-v1.na33.1.mm9.transcript.csv'
][organism]

result_file = [
    wd2 + 'human_transcriptomics_annot.tab',
    wd2 + 'mouse_transcriptomics_annot.tab'
][organism]

result_file_rma = [
    wd2 + 'human_transcriptomics_annot_rma.tab',
    wd2 + 'mouse_transcriptomics_annot_rma.tab'
][organism]

array_details_df = read_csv(array_details, comment='#', index_col='transcript_cluster_id')

# Read limma preprocessed data
array = read_csv([
    wd2 + 'human_transcriptomics.tab',
    wd2 + 'mouse_transcriptomics.tab'
][organism], sep='\t')

array['gene'] = [str(i).split(' // ')[1] if len(str(i).split(' // ')) > 1 else np.NaN for i in array_details_df.ix[array.index, 'gene_assignment']]
array = array.groupby('gene').min()

array.to_csv(result_file, sep='\t')

# Read rma normalised transcriptomics files
array = read_csv([
    wd + 'human_transcriptomics/human_transcriptomics_rma.csv',
    wd + 'mouse_transcriptomics/mouse_transcriptomics_rma.csv'
][organism], sep='\t', index_col=0)

array['gene'] = [str(i).split(' // ')[1] if len(str(i).split(' // ')) > 1 else np.NaN for i in array_details_df.ix[array.index, 'gene_assignment']]
array = array.groupby('gene').max()

array.to_csv(result_file_rma, sep='\t')