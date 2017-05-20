#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import rpathways
from scipy.stats import wilcoxon
from scipy_sugar.stats import quantile_gaussianize
from hlrcc.enrichment.gsea import gsea
from framed import load_cbmodel
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
fluxes = read_csv('./data/pfba_atp.csv', index_col=0).replace(np.nan, 0)
fluxes['delta'] = fluxes['UOK262'] - fluxes['UOK262pFH']


# -- Plotting
pal = dict(zip(*(['UOK262', 'UOK262pFH'], sns.light_palette('#34495e', 3, reverse=True).as_hex()[:-1])))

# ATP and biomass yield
plot_df = DataFrame([
    {
        'condition': 'UOK262pFH',
        'atp': abs(fluxes.ix['R_ATPM', 'UOK262pFH'] / fluxes.ix['R_EX_glc_e', 'UOK262pFH']),
        'atp_total': fluxes.ix['R_ATPM', 'UOK262pFH'],
        'biomass': abs(fluxes.ix['R_biomass_reaction', 'UOK262pFH'] / fluxes.ix['R_EX_glc_e', 'UOK262pFH'])
    },
    {
        'condition': 'UOK262',
        'atp': abs(fluxes.ix['R_ATPM', 'UOK262'] / fluxes.ix['R_EX_glc_e', 'UOK262']),
        'atp_total': fluxes.ix['R_ATPM', 'UOK262'],
        'biomass': abs(fluxes.ix['R_biomass_reaction', 'UOK262'] / fluxes.ix['R_EX_glc_e', 'UOK262'])
    }
])

order = ['UOK262', 'UOK262pFH']

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
g = sns.factorplot('condition', 'atp_total', data=plot_df, palette=pal, kind='bar', lw=0, order=order)
g.despine(right=False, top=False)
plt.gcf().set_size_inches(1, 2)
plt.ylabel('ATP production')
plt.xlabel('')
plt.savefig('./reports/atp_yield_bar.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'

sns.set(style='ticks', context='paper', font_scale=0.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in'})
g = sns.factorplot('condition', 'biomass', data=plot_df, palette=pal, kind='bar', lw=0, order=order)
g.despine(right=False, top=False)
plt.gcf().set_size_inches(1, 2)
plt.ylabel('Biomass yield (per mol Glucose)')
plt.xlabel('')
plt.savefig('./reports/biomass_yield_bar.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Plotting: Fluxes heatmaps
cmap = sns.diverging_palette(10, 220, sep=5, n=20, as_cmap=True)

# Flux pathway heatmap
plot_df = fluxes.ix[rpathways['glycolysis'] + rpathways['mitochondria']].drop(['R_CO2tm', 'R_ASPTAm', 'R_ASPTA'])
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

pal = dict(zip(*(['glycolysis', 'mitochondria'], sns.color_palette('Set2', n_colors=6).as_hex())))
rcol = Series({i: pal[p] for i in plot_df.index for p in rpathways if i in rpathways[p]}, name='')

plot_df.index = [i[2:] for i in plot_df.index]

sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.heatmap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, square=True)
plt.gcf().set_size_inches(.5, 4.5)
plt.setp(g.get_yticklabels(), rotation=0)
plt.setp(g.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# Exchange reactions
plot_df = fluxes.ix[{r for r in fluxes.index if r.startswith('R_EX_')}].sort('UOK262').dropna()
plot_df = plot_df[plot_df['delta'].abs() > 1e-5]

sns.set(style='white', context='paper', font_scale=.5, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.clustermap(plot_df, linewidths=.5, mask=(plot_df == 0), cmap=cmap, row_cluster=False, col_cluster=False)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
plt.gcf().set_size_inches(.5, 4)
plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.savefig('./reports/flux_heatmap_exchange.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot exported'


# -- Fluxomics enrichment
p_reactions = DataFrame([{'r': r, 'p': model.reactions[r].metadata['SUBSYSTEM']} for r in model.reactions if 'SUBSYSTEM' in model.reactions[r].metadata and model.reactions[r].metadata['SUBSYSTEM'] != ''])
p_reactions = p_reactions.groupby('p')['r'].agg(lambda x: set(x)).to_dict()

# Enrichment pathway analysis - Fluxomics
dataset = fluxes['delta'].replace(0, np.nan).dropna()
dataset = dict(zip(*(dataset.index, quantile_gaussianize(dataset))))

df_enrichment = [(c, len(sig.intersection(dataset)), gsea(dataset, sig, 1000)) for c, sig in p_reactions.items()]
df_enrichment = DataFrame([{'pathway': c, 'length': l, 'escore': es, 'pvalue': pval} for c, l, (es, pval) in df_enrichment]).dropna().set_index('pathway')
df_enrichment.sort(['escore']).to_csv('./files/fluxomics_pathways.csv', index=False)
print df_enrichment[df_enrichment['length'] > 2].sort('escore')
