#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import rpathways
from scipy.stats import wilcoxon
from pymist.enrichment.gsea import gsea
from framed import load_cbmodel
from protein_attenuation.utils import gkn
from pandas import read_csv, DataFrame, Series
from statsmodels.stats.multitest import multipletests


# -- Gene map
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])


# -- Import metabolic model
model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))


# -- Imports
conditions = ['UOK262', 'UOK262pFH']

ko_sampling, wt_sampling = [read_csv('./data/%s_sampling.txt' % c, sep='\t') for c in ['UOK262', 'UOK262pFH']]
fluxes = DataFrame(
    {r: {
        'UOK262': ko_sampling[r].median(),
        'UOK262pFH': wt_sampling[r].median(),
        'wilcoxon': wilcoxon(ko_sampling[r], wt_sampling[r])[1]
    } for r in ko_sampling}
).T
fluxes['fdr'] = multipletests(fluxes['wilcoxon'], method='fdr_bh')[1]
fluxes['delta'] = fluxes['UOK262'].abs() - fluxes['UOK262pFH'].abs()
print fluxes.sort(['delta'])


# -- Sampling plot
reactions = ['R_FUMm', 'R_MDHm', 'R_LEUTA', 'R_PGK', 'R_VALTA', 'R_GAPD', 'R_ADK1', 'R_ITCOALm']

ko_s = ko_sampling[reactions].unstack().reset_index()
ko_s.columns = ['Reaction', 'Index', 'Flux']
ko_s['Sample'] = 'KO'

wt_s = wt_sampling[reactions].unstack().reset_index()
wt_s.columns = ['Reaction', 'Index', 'Flux']
wt_s['Sample'] = 'WT'

plot_df = ko_s.append(wt_s)
pal = dict(zip(*(['KO', 'WT'], sns.light_palette('#34495e', 3).as_hex()[1:])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.FacetGrid(plot_df, col='Reaction', col_wrap=4, legend_out=True, sharey=False, aspect=1, size=1.2)
g.map(sns.boxplot, 'Sample', 'Flux', linewidth=.3, notch=True, fliersize=1, palette=pal, order=['KO', 'WT'])
g.map(plt.axhline, y=0, ls='-', lw=.3, alpha=.7, color='black')
g.set_titles('{col_name}')
g.set_ylabels('Flux rate (mmol/gDW/h)')
g.despine(trim=True)
plt.savefig('./reports/recon_sampling_plot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plotting: Fluxes heatmaps
cmap = sns.diverging_palette(10, 220, sep=5, n=20, as_cmap=True)

# Flux pathway heatmap
plot_df = fluxes.ix[rpathways['glycolysis'] + rpathways['mitochondria']].dropna()
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

pal = dict(zip(*(['glycolysis', 'mitochondria'], sns.color_palette('Set2', n_colors=6).as_hex())))
rcol = Series({i: pal[p] for i in plot_df.index for p in rpathways if i in rpathways[p]}, name='')

sns.set(style='white', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False, row_colors=rcol)
plt.gcf().set_size_inches(.5, 5)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# Heme flux pathway
plot_df = fluxes.ix[rpathways['heme']].dropna()
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

pal = dict(zip(*(['heme'], sns.color_palette('Set2', n_colors=6).as_hex())))
rcol = Series({i: pal[p] for i in plot_df.index for p in rpathways if i in rpathways[p]}, name='')

sns.set(style='white', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False, row_colors=rcol)
plt.gcf().set_size_inches(.5, 3)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap_heme.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# Exchange reactions
plot_df = fluxes.ix[{r for r in fluxes.index if r.startswith('R_EX_')}].sort('UOK262').dropna()
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False)
plt.gcf().set_size_inches(.5, 4)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap_exchange.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# -- Fluxomics enrichment
p_reactions = DataFrame([{'r': r, 'p': model.reactions[r].metadata['SUBSYSTEM']} for r in model.reactions if 'SUBSYSTEM' in model.reactions[r].metadata and model.reactions[r].metadata['SUBSYSTEM'] != ''])
p_reactions = p_reactions.groupby('p')['r'].agg(lambda x: set(x)).to_dict()

# Enrichment pathway analysis - Fluxomics
dataset = gkn(fluxes['delta'].replace(0, np.nan).dropna()).to_dict()

df_enrichment = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, 1000)) for c, sig in p_reactions.items()]
df_enrichment = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment]).dropna().set_index('pathway')
df_enrichment.sort(['escore']).to_csv('./files/fluxomics_pathways.csv', index=False)
print df_enrichment[df_enrichment['length'] > 2].sort('escore')
