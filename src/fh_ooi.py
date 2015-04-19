__author__ = 'emanuel'

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv

data_wd = '/Users/emanuel/Projects/data/'

samp = read_csv(data_wd + 'GSE26574/samplesheet.tab', sep='\t', index_col='Series_sample_id')
data = read_csv(data_wd + 'GSE26574/GSE26574_LOG2_RMA_HG-U133-Plus-2.0_ANNOT.tab', sep='\t', index_col='GeneSymbol')

data = data[samp.loc[[i in ['FH', 'FHN'] for i in samp['Phenotypes']], 'FileName']]

conditions = {
    'tumour': samp.loc[[i in ['FH'] for i in samp['Phenotypes']], 'FileName'].values,
    'normal': samp.loc[[i in ['FHN'] for i in samp['Phenotypes']], 'FileName'].values
}

genes = ['FH', 'CDH1', 'VIM', 'ZEB1', 'ZEB2']

# Barplot
data_df = DataFrame([(g, c, data.loc[g, s].median()) for g in genes for c, s in conditions.items()], columns=['gene', 'condition', 'expression'])
sns.set_style('white')
g = sns.factorplot('gene', 'expression', 'condition', data_df, kind='bar', palette=sns.color_palette('Paired'), aspect=1.25, legend=False)
sns.despine()
g.despine(offset=10, trim=True)
g.set_axis_labels('', 'Expression (log2)')
plt.legend(loc='upper right')
plt.savefig('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/plots/' + 'ooi_et_al_expression.pdf')
plt.close('all')

# Boxplot
data_df = DataFrame([(g, c, v) for g in genes for c, s in conditions.items() for v in data.loc[g, s].values], columns=['gene', 'condition', 'expression'])

sns.set_style('white')
g = sns.factorplot('gene', 'expression', 'condition', data_df, kind='box', palette=sns.color_palette('Paired'), aspect=1.25, legend=False)
sns.despine()
g.despine(offset=10, trim=True)
g.set_axis_labels('', 'Expression (log2)')
plt.legend(loc='upper right')
plt.savefig('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/plots/' + 'ooi_et_al_expression_boxplot.pdf')
plt.close('all')

sample_data = data.ix[genes]
sample_data.columns = ['normal' if c in conditions['normal'] else 'tumour' for c in data.columns]
sample_data.to_csv('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/plots/' + 'ooi_et_al_expression.csv')