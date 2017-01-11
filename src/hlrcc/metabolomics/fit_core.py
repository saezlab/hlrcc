#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
from scipy.stats import spearmanr, pearsonr
from pymist.reader.sbml_reader import read_sbml_model
from pandas import read_csv, DataFrame, Series
from pymist.simulation import min_differences


# -- Imports
conditions = ['UOK262', 'UOK262pFH']

# Import metabolite map
m_map = read_csv('./files/metabolites_map.txt', sep='\t', index_col=0)
m_map = m_map.to_dict()['metabolite']

# Metabolic model
model = read_sbml_model('./files/recon1.xml')
model.remove_b_metabolites()

# CORE metabolomics (mmol/gDW/h)
core = read_csv('./data/uok262_metabolomics_core_processed.txt', sep='\t', index_col=0).dropna().set_index('exchange')
core = core[core['fdr'] < .05]

# Medium
medium = read_csv('./files/DMEM_41966_medium.txt', sep='\t').dropna()

# O2 consumption (mmol/gDW/h)
o2_exch = 'R_EX_o2_e_'

# -- Estimate exchange reactions
# Set condition
predictions = {}
for condition in conditions:
    # condition = conditions[0]
    print '[INFO] Condition: ', condition

    # Get condition exchange rates
    exchange_rates = core[condition].dropna()

    # -- Fit medium
    for metabolite in exchange_rates.index:
        # metabolite = exchange_rates.index[0]
        # Reactions not used
        for reaction in ['R_EX_retpalm', 'R_EX_retpalm_e_']:
            model.set_constraint(reaction, 0)

        # Deepcopy metabolic model
        reduced_model = model.deepcopy()

        # Add FH KO to tumour cell lines
        if condition == 'UOK262':
            reduced_model.set_constraint('R_FUM', 0, 0)
            reduced_model.set_constraint('R_FUMm', 0, 0)
            print '[INFO] R_FUM blocked'

        # Block exhange fluxes not present in the medium or CORE experiment
        for reaction in reduced_model.get_exchanges(check_matrix=True):
            if (reaction not in exchange_rates.drop(metabolite).index) and (reaction != o2_exch) and (reaction not in medium['exchange'].values):
                reduced_model.set_constraint(reaction, lower_bound=0)
        print '[INFO] Exchange reactions constrained'

        # Minimise differences between simulations and measurements
        constrains = {}
        for r in reduced_model.get_exchanges(check_matrix=True):
            if r in exchange_rates.drop(metabolite).index:
                constrains[r] = exchange_rates.ix[r, condition]
            # elif r == 'R_EX_o2_e_':
            #     constrains['R_EX_o2_e_'] = core_o2[condition]

        fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model, check_matrix=True)

        # Constraint metabolic model lower bound
        for reaction in fitted_medium:
            if reaction in exchange_rates.drop(metabolite).index:
                reduced_model.set_constraint(reaction, lower_bound=fitted_medium[reaction], upper_bound=fitted_medium[reaction])

        # Run pFBA
        internal_reactions = set(reduced_model.reactions) - set(reduced_model.get_exchanges(check_matrix=True))
        constrains = {r: 0 for r in internal_reactions}
        fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model)

        # Save fitted medium measurements
        predictions['%s %s' % (condition, metabolite)] = fitted_medium[metabolite] if metabolite in fitted_medium else 0
        print metabolite, predictions['%s %s' % (condition, metabolite)]

predictions = DataFrame(Series(predictions), columns=['Predicted'])

predictions['Measured'] = [core.ix[i.split(' ')[1], i.split(' ')[0]] for i in predictions.index]
print predictions


# -- Evaluate predictions
sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.jointplot(
    'Measured', 'Predicted', predictions, 'reg', color='#34495e', joint_kws={'ci': None, 'scatter_kws': {'s': 40, 'edgecolor': 'w', 'linewidth': .5}},
    marginal_kws={'hist': False, 'rug': True}, annot_kws={'template': 'Spearman: {val:.2g}, p-value: {p:.1e}'}, space=0,
    stat_func=spearmanr
)
plt.axhline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
plt.axvline(0, ls='-', lw=0.3, c='#95a5a6', alpha=.5)
g.plot_marginals(sns.kdeplot, shade=True, color='#34495e')
g.set_axis_labels('Measured (umol/ugDW/h)', 'Predicted (umol/ugDW/h)')
plt.gcf().set_size_inches(3, 3)
plt.savefig('./reports/fitted_medium_loo_cor.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'


# -- Fit metabolic model with all CORE rates
mediums = {}
for condition in conditions:
    # condition = conditions[0]
    print '[INFO] Condition: ', condition

    # Get condition exchange rates
    exchange_rates = core[condition].dropna()

    # -- Fit medium
    # Reactions not used
    for reaction in ['R_EX_retpalm', 'R_EX_retpalm_e_']:
        model.set_constraint(reaction, 0)

    # Deepcopy metabolic model
    reduced_model = model.deepcopy()

    # Add FH KO to tumour cell lines
    if condition == 'UOK262':
        reduced_model.set_constraint('R_FUM', 0, 0)
        reduced_model.set_constraint('R_FUMm', 0, 0)
        print '[INFO] R_FUM blocked'

    # Block exhange fluxes not present in the medium or CORE experiment
    for reaction in reduced_model.get_exchanges(check_matrix=True):
        if (reaction not in exchange_rates.index) and (reaction != o2_exch) and (reaction not in medium['exchange'].values):
            reduced_model.set_constraint(reaction, lower_bound=0)
    print '[INFO] Exchange reactions constrained'

    # Minimise differences between simulations and measurements
    constrains = {}
    for r in reduced_model.get_exchanges(check_matrix=True):
        if r in exchange_rates.index:
            constrains[r] = exchange_rates.ix[r, condition]

    fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model, check_matrix=True)

    # Constraint metabolic model lower bound
    for reaction in fitted_medium:
        if reaction in exchange_rates.index:
            reduced_model.set_constraint(reaction, lower_bound=fitted_medium[reaction], upper_bound=fitted_medium[reaction])

    # Run pFBA
    internal_reactions = set(reduced_model.reactions) - set(reduced_model.get_exchanges(check_matrix=True))
    constrains = {r: 0 for r in internal_reactions}
    fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model)
    
    # Save fitted medium measurements
    mediums[condition] = fitted_medium

mediums = DataFrame(mediums)
print mediums

# Plot
plot_df = mediums.ix[mediums.eval('-'.join(conditions)).abs().sort_values(inplace=False, ascending=False).index].unstack().reset_index().dropna()
plot_df.columns = ['condition', 'exchange', 'rate']
plot_df['metabolite'] = [model.metabolites[model.get_reactants(e)[0]] for e in plot_df['exchange']]

pallete = sns.light_palette('#34495e', 3)[1:]

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
g = sns.factorplot('rate', 'metabolite', data=plot_df, hue='condition', palette=pallete, legend=True, legend_out=True, aspect=.5, size=3, scale=.5)
plt.axvline(0, c='#95a5a6', lw=.3, alpha=.7, ls='-')
plt.xlabel('umol/ugDW/h')
plt.ylabel('')
sns.despine(trim=True)
plt.savefig('./reports/fitted_medium.pdf', bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
