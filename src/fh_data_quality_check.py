__author__ = 'emanuel'

import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
from pandas import read_csv, Series, DataFrame
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


def run_cv(dataset, name, title, conditions, plot_dir):
    for (i, c) in conditions:
        plt.subplot(1, 2, i+1)
        dataset_cv = dataset[c].std(axis=1) / dataset[c].mean(axis=1)
        plt.hist(dataset_cv.dropna(), 20, histtype='stepfilled', alpha=.7)
        sns.despine()
        plt.title(title + ': ' + c + ' - Coeficient variation')

    save_plot(plot_dir + name + 'coefficient_variation.pdf')


def run_pca(dataset, name, title, plot_dir):
    pca = PCA(n_components=2)
    dataset_pca = DataFrame(pca.fit_transform(dataset.dropna().T), index=dataset.columns)

    [plt.scatter(dataset_pca.loc[c, 0], dataset_pca.loc[c, 1], c=col_pal_1[i], label=c) for (i, c) in conditions]
    sns.despine()

    plt.xlabel('PC 1 (%.2f)' % pca.explained_variance_ratio_[0])
    plt.ylabel('PC 2 (%.2f)' % pca.explained_variance_ratio_[1])
    plt.legend(frameon=True, fancybox=True)

    plt.title(title)

    save_plot(plot_dir + name + 'pca.pdf', 5, 5)


def save_plot(path, w=12., h=7.):
    fig = plt.gcf()
    fig.set_size_inches(w, h)
    fig.savefig(path)
    plt.close()
    print '[INFO] Plot saved: ' + path


# Configure
data_wd, plot_dir = '/Users/emanuel/Projects/data/fh_cells/', '/Users/emanuel/Projects/projects/pymist/reports/fh_cells/'

# Colour pallete
col_pal_1 = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
conditions = [(0, 'fh_ko'), (1, 'fh_wt')]

# Import data-sets
ss = read_csv(data_wd + 'fh_samplesheet.tab', sep='\t', index_col=0)

pp_human_all = read_csv(data_wd + 'human_phosphoproteomics/b1368p100_phosho_human.tab', sep='\t', index_col=0)
tp_human_all = read_csv(data_wd + 'human_proteomics/b1368p100_protein_human.tab', sep='\t', index_col=0)

# Assemble reduced data-set
pp_human = pp_human_all[ss.loc[(ss['organism'] == 'human') & (ss['type'] == 'pp')].index]
pp_human.columns = ss.loc[(ss['organism'] == 'human') & (ss['type'] == 'pp')]['condition']

tp_human = tp_human_all[ss.loc[(ss['organism'] == 'human') & (ss['type'] == 'tp')].index]
tp_human.columns = ss.loc[(ss['organism'] == 'human') & (ss['type'] == 'tp')]['condition']

# Plot original distributions
pp_human.hist(); plt.show();
tp_human.hist(); plt.show();

# Log2 transform
pp_human = np.log2(pp_human)
pp_human[np.isinf(pp_human)] = np.NaN

tp_human = np.log2(tp_human)
tp_human[np.isinf(tp_human)] = np.NaN

# Plot log2 intensities
sns.set(style='white')
sns.violinplot(pp_human.dropna(), vert=True, pallete='Set2')
sns.despine()
plt.title('Phosphoproteomics')
plt.xlabel('peptides intensities (log2)')
save_plot(plot_dir + 'human_phospho_log2_replicatess_distribution_violin_plots.pdf', 10, 5)

sns.set(style='white')
sns.violinplot(tp_human.dropna(), vert=True, pallete='Set2')
sns.despine()
plt.title('Proteomics')
plt.xlabel('peptides intensities (log2)')
save_plot(plot_dir + 'human_proteomics_log2_replicatess_distribution_violin_plots.pdf', 5, 5)

# NAs counts
pp_nan_count = np.isnan(pp_human).sum(axis=1)
tp_nan_count = np.isnan(tp_human).sum(axis=1)

plt.hist(pp_nan_count, 8, histtype='stepfilled', alpha=.7)
sns.despine()
plt.title('Phosphoproteomics')
plt.xlabel('NaN counts')
save_plot(plot_dir + 'human_phospho_nan_counts.pdf', 5, 5)

plt.hist(tp_nan_count, 8, histtype='stepfilled', alpha=.7)
sns.despine()
plt.title('Proteomics')
plt.xlabel('NaN counts')
save_plot(plot_dir + 'human_proteomics_nan_counts.pdf', 5, 5)

# Coefficient variation
run_cv(pp_human, 'human_phospho_', 'phosphoproteomics', conditions, plot_dir)
run_cv(tp_human, 'human_protein_', 'proteomics', conditions, plot_dir)

# PCA
run_pca(pp_human, 'human_phospho_', 'phosphoproteomics', plot_dir)
run_pca(tp_human, 'human_protein_', 'proteomics', plot_dir)

# Dendogram
distxy = squareform(pdist(pp_human.dropna().T, metric='euclidean'))
dendog = dendrogram(linkage(distxy, method='complete'))
plt.show()

# Correlation matrix
f, ax = plt.subplots(figsize=(9, 9))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.corrplot(pp_human, cmap=cmap, ax=ax)
f.tight_layout()
plt.title('Phosphoproteomics')
save_plot(plot_dir + 'human_phosphop_coplot.pdf', 12, 12)

f, ax = plt.subplots(figsize=(9, 9))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
sns.corrplot(tp_human, cmap=cmap, ax=ax)
f.tight_layout()
plt.title('Proteomics')
save_plot(plot_dir + 'human_proteomics_coplot.pdf', 12, 12)