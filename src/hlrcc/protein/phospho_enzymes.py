#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from framed import load_cbmodel
from more_itertools import unique_everseen
from scipy.stats.stats import spearmanr, pearsonr
from pandas import read_csv, Series, DataFrame

# -- Imports
conditions = ['KO', 'WT']
conditions_map = {'UOK262': 'KO', 'UOK262pFH': 'WT'}

# Proteomics
proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')
phosphoproteomics = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')

# Fluxomics
fluxes = read_csv('./data/pfba_atp.csv', index_col=0).replace(np.nan, 0)
fluxes['delta'] = fluxes['UOK262'] - fluxes['UOK262pFH']
fluxes.columns = [conditions_map[c] if c in conditions_map else c for c in fluxes]

# Regulatory p-sites
r_sites = read_csv('./files/Regulatory_sites.txt', sep='\t')
r_sites = r_sites[r_sites['ORGANISM'] == 'human']
r_sites = r_sites[[i.endswith('-p') for i in r_sites['MOD_RSD']]]
r_sites['res'] = ['%s_%s' % (g, p.split('-')[0]) for g, p in r_sites[['GENE', 'MOD_RSD']].values]


# -- Import metabolic model
gmap = read_csv('./files/non_alt_loci_set.txt', sep='\t')
gmap['hgsn'] = ['G_' + i.replace(':', '_') for i in gmap['hgnc_id']]
gmap = gmap.groupby('hgsn')['symbol'].agg(lambda x: list(x)[0])

model = load_cbmodel('./files/recon2.2.xml', flavor='cobra')
model.detect_biomass_reaction()
model.remove_metabolite('M_biomass_c')
model.add_reaction_from_str('R_ATPM: M_h2o_c + M_atp_c --> M_adp_c + M_pi_c + M_h_c')

g_enzymes = {gmap[g] for g in model.genes if g in gmap}

g_to_r = model.gene_to_reaction_lookup()
g_to_r = {gmap[g]: g_to_r[g] for g in g_to_r if g in gmap}


# --
plot_df = DataFrame({p: {
    'proteomics': proteomics[p.split('_')[0]],
    'phosphoproteomics': phosphoproteomics[p]
} for p in phosphoproteomics.index if p.split('_')[0] in proteomics}).T
plot_df['protein'] = [p.split('_')[0] for p in plot_df.index]
plot_df = plot_df[[i in g_enzymes for i in plot_df['protein']]]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'proteomics', 'phosphoproteomics', plot_df, 'reg', color='#34495e', space=0,
    joint_kws={'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5, 'alpha': .5}, 'fit_reg': False},
    marginal_kws={'hist': False, 'rug': False}, stat_func=spearmanr,
    annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}', 'loc': 4}
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Proteomics (log2 FC)', 'Phosphoproteomics (log2 FC)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/proteomics_phospho_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'


# -- Regress-out protein abundance
lm = sm.OLS(plot_df['phosphoproteomics'], sm.add_constant(plot_df['proteomics'])).fit()
residuals = lm.resid.sort_values()
print lm.summary()

# Plot
plot_df = residuals.reset_index()
plot_df.columns = ['psite', 'change']
plot_df['abs_change'] = plot_df['change'].abs()
plot_df = plot_df.sort_values('abs_change', ascending=False).head(20)
plot_df = plot_df.sort_values('change', ascending=False)

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.stripplot(y='psite', x='change', data=plot_df, color='#bfbfbf', s=4)
plt.axvline(0, ls='-', lw=0.3, c='black', alpha=.5)
plt.xlabel('Fold-change (protein normalised)')
plt.ylabel('')
plt.legend(loc=0, title='Type')
# plt.xticks(np.arange(-1.5, 2.5, .5))
plt.gcf().set_size_inches(1.5, 2.5)
sns.despine(trim=True)
plt.title('Metabolic enzymes\nphosphorylation-sites')
plt.savefig('./reports/psite_regression_stripplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# --
plot_df = DataFrame([{
    'Psite': p,
    'Enzyme': p.split('_')[0],
    'Reaction': r,
    'phospho (logfc)': phosphoproteomics[p],
    'flux (delta)': fluxes.ix[r, 'delta']
 } for p in residuals.index for r in g_to_r[p.split('_')[0]] if r in fluxes.index]).sort_values('flux (delta)')
plot_df = plot_df[plot_df['Reaction'] != 'R_ATPS4m']
plot_df = plot_df[plot_df['flux (delta)'].abs() > 1e-1]
plot_df['flux (KO)'] = [fluxes.ix[i, 'KO'] for i in plot_df['Reaction']]
plot_df['flux (WT)'] = [fluxes.ix[i, 'WT'] for i in plot_df['Reaction']]
plot_df.to_csv('./data/fluxes_phospho_table.csv', index=False)
print plot_df.sort_values('flux (delta)')

plot_df[[i in set(r_sites['res']) for i in plot_df['Psite']]]
r_sites.loc[[i in set(plot_df['Psite']) for i in r_sites['res']], ['res', 'NOTES']]

print r_sites.ix[1492]
print r_sites.ix[1492, 'NOTES']

order = list(unique_everseen(plot_df['Enzyme']))
pal = dict(zip(*(order, sns.diverging_palette(10, 220, sep=5, n=len(order)+1).as_hex()[:-1])))
# pal = dict(zip(*(order, sns.light_palette('#34495e', len(order)+1, reverse=True).as_hex()[:-1])))

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.lmplot(
    'flux (delta)', 'phospho (logfc)', plot_df, 'Enzyme', fit_reg=False, palette=pal, size=3,
    scatter_kws={'s': 50, 'edgecolor': 'w', 'linewidth': .5}, legend_out=False, legend=False
)
plt.axhline(0, ls='-', lw=.3, c='gray')
plt.axvline(0, ls='-', lw=.3, c='gray')
plt.xlabel('Fluxomics (mean difference)')
plt.ylabel('Phosphorylation-site (log2 FC)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Enzymes')
plt.gcf().set_size_inches(2.5, 2.5)
plt.savefig('./reports/fluxomics_phospho_jointplot.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Corr plotted!'
