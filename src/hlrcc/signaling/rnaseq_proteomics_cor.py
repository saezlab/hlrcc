import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from pandas import DataFrame, Series, read_csv
from pymist.utils.map_peptide_sequence import read_uniprot_genename


# -- Imports
# Ensembl ids
ensembl = read_csv('%s/files/ensembl_id_map_human.txt' % wd, sep='\t', index_col=0).dropna().to_dict()['HGNC symbol']

# Uniprot ids
human_uniprot = read_uniprot_genename()
print '[INFO] Uniprot human protein: ', len(human_uniprot)

# Rna-seq
rnaseq = read_csv('%s/data/UOK262_rnaseq_logfc.txt' % wd, sep='\t', index_col=0)
rnaseq['name'] = [ensembl[i] if i in ensembl else np.NaN for i in rnaseq.index]
rnaseq = rnaseq.dropna()
rnaseq = rnaseq.set_index('name')

# Proteomics
proteomics = read_csv('%s/data/uok262_proteomics_logfc.txt' % wd, sep='\t', index_col=0)
proteomics['name'] = [human_uniprot[i.split('_')[0]][0] if i.split('_')[0] in human_uniprot else np.nan for i in proteomics.index]
proteomics = proteomics.dropna()
proteomics = proteomics.set_index('name')


# -- Correlate
ov_meas = set(rnaseq.index).intersection(proteomics.index)

plot_df = DataFrame({'rnaseq': rnaseq.ix[ov_meas, 'logFC'], 'proteomics': proteomics.ix[ov_meas, 'logFC']})

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'rnaseq', 'proteomics', plot_df, 'reg', color='#34495e', joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}},
    marginal_kws={'hist': False, 'rug': False}, annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4}, ylim=[-5, 5], space=0
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Transcriptomics (log2 FC)', 'Proteomics (log2 FC)')
plt.savefig('%s/reports/rnaseq_proteomics_jointplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
